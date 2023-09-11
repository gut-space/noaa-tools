# This script processes images received from NOAA satellites

from noaatools.constants import DEG2RAD, RAD2DEG, Ellipsoid, Method, NOAA_PROCESSING_DELAY, RE, AVHRR_FOV, ellipsoid_wgs84
from noaatools import export_js
import sys
from datetime import datetime, timezone, timedelta
from math import atan, atan2, sqrt, pi, sin, cos, asin, acos, tan
from typing import Tuple

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.api import jday, Satrec

import numpy as np
from pymap3d import ecef

sys.path.append('.')

# Nice conversions: https://github.com/skyfielders/python-skyfield/blob/master/skyfield/sgp4lib.py
# Good explanation: https://stackoverflow.com/questions/8233401/how-do-i-convert-eci-coordinates-to-longitude-latitude-and-altitude-to-display-o


def julianDateToGMST(jd, fr):
    """
    Converts Julian date (expressed at two floats) to GMST (Greenwich Mean Sidereal Time).

    Parameters:
    jd : float - Julian date full integer + 0.5
    fr : float - fractional part of the Julian date

    Returns
    =======
    A single floating point representing a GMST, expressed in degrees (0...359.99999).

    This calculation takes into consideration the precession, but not nutation.

    Source: https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#GMST
    """
    T0 = 2451545.0  # J2000, 2000-Jan-01 12h UT1 as Julian date

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - T0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    T = d / 36525.0

    # Calculate GMST (in seconds at UT1=0)
    gmst = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 0.0000062 * T * T * T

    # Let's truncate this and return the value in degrees.
    # This is clearly broken.
    return (gmst % 24) * (15 / 3600.0)


def julianDateToGMST2(jd: float, fr: float) -> Tuple[float, float]:
    """
    Converts Julian date (expressed at two floats) to GMST (Greenwich Mean Sidereal Time 1982).

    Parameters:
    jd : float - Julian date full integer + 0.5
    fr : float - fractional part of the Julian date

    Returns
    =======
    A tuple with two values:
    theta - single floating point representing a GMST, expressed in radians
    theta_dot - unknown

    This calculation takes into consideration the precession, but not nutation.

    Source: https://github.com/skyfielders/python-skyfield/blob/master/skyfield/sgp4lib.py
    - theta_GSMT1982 function

    This angle defines the difference between the idiosyncratic True
    Equator Mean Equinox (TEME) frame of reference used by SGP4 and the
    more standard Pseudo Earth Fixed (PEF) frame of reference.

    From AIAA 2006-6753 Appendix C.
    """
    tau = 6.283185307179586476925287

    _second = 1.0 / (24.0 * 60.0 * 60.0)

    T0 = 2451545.0  # J2000, 2000-Jan-01 12h UT1 as Julian date

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - T0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    t = d / 36525.0

    # Don't undersran
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = ((jd + fr) % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau
    return theta, theta_dot


def longitude_trunc(lon: float) -> float:
    """
    Makes sure the longitude is within <-pi ... pi> range.

    Parameters
    ==========
    lon - longitude expressed in radians (may be any value)

    Returns
    =======
    normalized longitude in <-pi..pi> range. Note that both -pi and pi are accepted.
    """

    if (lon <= pi) and (lon >= -pi):
        # Don't do any conversion if it's not necessary. Avoid conversion errors if necessary
        return lon
    return (lon + pi) % (2 * pi) - pi


def teme2geodetic_spherical(x: float, y: float, z: float, t: datetime):
    """
    Converts TEME/ECI coords (x,y,z - expressed in km) to LLA (longitude, lattitude, altitude).
    This function assumes the Earth is completely round.

    The calculations here are based on T.S. Kelso's excellent paper "Orbital Coordinate Systems, Part III
    https://celestrak.com/columns/v02n03/.

    Parameters
    ==========
    x,y,z : float - coordates in TEME (True Equator Mean Equinoex) version of ECI (Earth Centered Intertial) coords system.
            This is the system that's produced by SGP4 models.
    t : datetime

    Returns
    =======

    lat, lon, alt - latitude, longitude (degrees), altitude (km)
    """

    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    gmst = julianDateToGMST2(jd, fr)[0]

    lat = atan2(z, sqrt(x * x + y * y))  # phi
    lon = atan2(y, x) - gmst  # lambda-E
    lon = longitude_trunc(lon)
    alt = sqrt(x * x + y * y + z * z) - RE  # h

    # TODO: convert this to radians and use radians everywhere.
    return lat * RAD2DEG, lon * RAD2DEG, alt


def teme2geodetic_oblate(x: float, y: float, z: float, t: datetime, ellipsoid: Ellipsoid):
    """
    Converts TEME/ECI coords (x,y,z - expressed in km) to LLA (longitude, lattitude, altitude).
    ellipsoid is Earth ellipsoid to be used (e.g. ellipsoid_wgs84).

    The calculations here are based on T.S. Kelso's excellent paper "Orbital Coordinate Systems, Part III
    https://celestrak.com/columns/v02n03/.

    Parameters
    ==========
    x,y,z : float - coordates in TEME (True Equator Mean Equinoex) version of ECI (Earth Centered Intertial) coords system.
            This is the system that's produced by SGP4 models.
    t : datetime - time of the observation
    ellipsoid: Ellipsoid - an Earth exlipsoid specifying Earth oblateness, e.g. Ellipsoid_wgs84. Two params are used from it:
            a and inverse of f. Both must be specified in kms

    Returns
    =======
    lat, lon, alt - latitude, longitude (both in degrees), alt (in km)
    """

    # First, we need to do some basic calculations for Earth oblateness
    a = ellipsoid.a
    f = 1.0 / ellipsoid.finv
    b = a * (1 - 1.0 / f)
    e2 = f * (2 - f)

    phii = 1  # This is the starting value for initial iteration

    # There should be a check on |phii - phi| value, but let's always do 5 iterations. Good enough for now.
    for _ in range(1, 5):

        C = 1 / (sqrt(1 - e2 * pow(sin(phii), 2)))
        # This is not explicitly stated on clestrak page, but it's shown on a diagram.
        R = sqrt(x * x + y * y)
        phi = atan2(z + a * C * e2 * sin(phii), R)
        h = R / (cos(phi)) - a * C

        phii = phi

    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    gmst = julianDateToGMST2(jd, fr)[0]

    lon = atan2(y, x) - gmst  # lambda-E
    lon = longitude_trunc(lon)

    return phi * RAD2DEG, lon * RAD2DEG, h


def teme2geodetic_pymap3d(x: float, y: float, z: float, t: datetime, ell=None):
    """
    Converts TEME/ECI coordinates to geodetic, using pymap3d library.
    For details, see https://github.com/geospace-code/pymap3d

    Parameters
    ==========
    x,y,z : float - coordates in TEME (True Equator Mean Equinoex) version of ECI (Earth Centered Intertial) coords system.
            This is the system that's produced by SGP4 models.
    t : datetime - time of the observation
    ell: unknown - haven't figured this out yet.

    Returns
    =======
    lat, lon, alt - latitude, longitude (in degrees), alt (in km)
    """

    # Short version - whole conversion in one go
    # lat, lon, alt = ecef.eci2geodetic(x*1000, y*1000, z*1000, t)
    # return lat, lon, alt/1000.0

    # print("teme[x,y,z]=%f, %f, %f" % (x, y, z))
    xecef, yecef, zecef = ecef.eci2ecef(np.array([x * 1000]), np.array([y * 1000]), np.array([z * 1000]), t)
    # print("ecef[x,y,z]=%f, %f, %f" % (xecef, yecef, zecef))

    # True = we want the response in degrees
    lat, lon, alt = ecef.ecef2geodetic(xecef, yecef, zecef, ell, True)
    # print("lla = %f, %f, %f" % (lat, lon, alt))
    return lat, lon, alt / 1000.0


def get_ssp(lla):
    return [lla[0], lla[1], 0]


def calc_azimuth(p1, p2):
    """ Calculates azimuth from point 1 to point 2.
    Point - an array 3 of floats (LLA)
    Returns azimuth in degrees

    Source: http://edwilliams.org/avform.htm#Crs
    """

    lat1 = p1[0] * DEG2RAD
    lon1 = -p1[1] * DEG2RAD
    lat2 = p2[0] * DEG2RAD
    lon2 = -p2[1] * DEG2RAD

    d = 2 * asin(sqrt((sin((lat1 - lat2) / 2))**2 + cos(lat1) * cos(lat2) * (sin((lon1 - lon2) / 2))**2))

    tc1 = acos((sin(lat2) - sin(lat1) * cos(d)) / (sin(d) * cos(lat1)))

    tc1 = atan2(sin(lon1 - lon2) * cos(lat2), cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2))
    if (tc1 < 0):
        tc1 += 2 * pi
    if (tc1 > 2 * pi):
        tc1 -= 2 * pi
    return tc1 * RAD2DEG


def calc_swath(alt, nu):
    """This calculates the swath width, given the altitude (alt, in km) of the sat and camera angle (nu, in radians).
        Returns swath in km"""

    # Convert to radians first.
    nur = nu * DEG2RAD

    # Ok, this is an overly simplified approximation. It neglects the Earth curvature.
    # return alt*tan(nu)

    # Source Wertz "Mission geometry", pg. 420.

    # rho is an angle between two lines: (sat - tangential to Earth) and (sat - Earth center)
    rho = asin(RE / (RE + alt))

    epsilon = acos(sin(nur) / sin(rho))
    lam = pi / 2 - nur - epsilon
    swath = RE * lam
    print("calc_swath(alt=%f nu=%f/%f) => rho=%f/%f epsilon=%f/%f, lambda= %f/%f => swath=%f [km]" %
          (alt, nu, nur, rho, rho * RAD2DEG, epsilon, epsilon * RAD2DEG, lam, lam * RAD2DEG, swath))

    return swath


def radial_distance(lat1, lon1, bearing, distance):
    """
    Return final coordinates (lat2,lon2) [in degrees] given initial coordinates
    (lat1,lon1) [in degrees] and a bearing [in degrees] and distance [in km]

    Based on this:
    https://stackoverflow.com/questions/877524/calculating-coordinates-given-a-bearing-and-a-distance
    """

    rlat1 = lat1 * DEG2RAD
    rlon1 = lon1 * DEG2RAD
    rdistance = distance / RE  # normalize linear distance to radian angle
    rbearing = bearing * DEG2RAD

    rlat = asin(sin(rlat1) * cos(rdistance) + cos(rlat1) * sin(rdistance) * cos(rbearing))

    if cos(rlat) == 0 or abs(cos(rlat)) < 0.00000001:  # Endpoint a pole
        rlon = rlon1
    else:
        rlon = ((rlon1 + asin(sin(rbearing) * sin(rdistance) / cos(rlat)) + pi) % (2 * pi)) - pi

    return (rlat * RAD2DEG, rlon * RAD2DEG)


def calc_distance(lat1, lon1, lat2, lon2):
    """
    Calculates distance between two (lat,lon) points. Return value is in km.
    """
    rlat1 = lat1 * DEG2RAD
    rlon1 = lon1 * DEG2RAD
    rlat2 = lat2 * DEG2RAD
    rlon2 = lon2 * DEG2RAD

    d = 2 * asin(sqrt((sin((rlat1 - rlat2) / 2))**2 + cos(rlat1) * cos(rlat2) * (sin((rlon1 - rlon2) / 2))**2))

    return d * RE


def azimuth_add(az, delta):
    """ Adds delta to specified azimuth. Does the modulo 360 arithmetic"""

    return (az + delta) % 360.0


def teme2geodetic(method: Method, x: float, y: float, z: float, t: datetime):
    if method == Method.SPHERICAL:
        return teme2geodetic_spherical(x, y, z, t)
    if method == Method.OBLATE:
        return teme2geodetic_oblate(x, y, z, t, ellipsoid_wgs84)
    if method == Method.PYMAP3D:
        return teme2geodetic_pymap3d(x, y, z, t)
    raise Exception("Invalid calculation method: %s" % method)


def georef(method: Method, tle1: str, tle2: str, aos_txt: str, los_txt: str):
    """ This is a naive georeferencing method:
        - calculates the sat location at AOS and LOS points (using )
    then calculates distance between them. """

    # Convert date as a string datetime. Make sure to use UTC rather than the default (local timezone)
    d1 = datetime.fromisoformat(aos_txt).replace(tzinfo=timezone.utc)
    d2 = datetime.fromisoformat(los_txt).replace(tzinfo=timezone.utc)

    print("AOS time: %s" % d1)
    print("LOS time: %s" % d2)

    # STEP 1: Calculate sat location at AOS and LOS

    # Old sgp4 API 1.x used this approach, which is not recommended anymore.
    # sat_old = twoline2rv(tle1, tle2, wgs72)
    # pos1_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    # pos2_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)

    # This approach uses new API 2.x which gives a slightly different results.
    # In case of NOAA, the position is off by less than milimeter
    sat = Satrec.twoline2rv(tle1, tle2)
    jd1, fr1 = jday(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    jd2, fr2 = jday(d2.year, d2.month, d2.day, d2.hour, d2.minute, d2.second)

    # Take sat processing/transmission delay into consideration. At AOS time the signal received
    # was already NOAA_PROCESSING_DELAY seconds old.
    fr1 -= NOAA_PROCESSING_DELAY / 86400.0
    fr2 -= NOAA_PROCESSING_DELAY / 86400.0

    _, pos1, _ = sat.sgp4(jd1, fr1)  # returns error, position and velocity - we care about position only
    _, pos2, _ = sat.sgp4(jd2, fr2)

    # Delta between a point and a point+delta (the second delta point is used to calculate azimuth)
    DELTA = 30.0

    _, pos1delta, _ = sat.sgp4(jd1, fr1 + DELTA / 86400.0)
    _, pos2delta, _ = sat.sgp4(jd2, fr2 + DELTA / 86400.0)

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

    # Ok, we have sat position at time of AOS and LOS returned by SGP4 models. The tricky part here is those are in
    # Earth-Centered Intertial (ECI) reference system. The model used is TEME (True equator mean equinox).

    # Convert AOS coordinates to LLA
    aos_lla = teme2geodetic(method, pos1[0], pos1[1], pos1[2], d1)

    # Now caluclate a position for AOS + 30s. Will use it to determine the azimuth
    d1delta = d1 + timedelta(seconds=30.0)
    aos_bis = teme2geodetic(method, pos1delta[0], pos1delta[1], pos1delta[2], d1delta)
    aos_az = calc_azimuth(aos_lla, aos_bis)

    print("AOS converted to LLA is lat=%f long=%f alt=%f, azimuth=%f" % (aos_lla[0], aos_lla[1], aos_lla[2], aos_az))

    # Now do the same for LOS
    los_lla = teme2geodetic(method, pos2[0], pos2[1], pos2[2], d2)

    # Let's use a point 30 seconds later. AOS and (AOS + 30s) will determine the azimuth
    d2delta = d2 + timedelta(seconds=30.0)
    los_bis = teme2geodetic(method, pos2delta[0], pos2delta[1], pos2delta[2], d2delta)
    los_az = calc_azimuth(los_lla, los_bis)

    print("LOS converted to LLA is lat=%f long=%f alt=%f azimuth=%f" % (los_lla[0], los_lla[1], los_lla[2], los_az))

    # STEP 3: Find image corners. Here's an algorithm proposal:
    #
    # 1. calculate satellite flight azimuth AZ
    #    https://en.wikipedia.org/wiki/Great-circle_navigation
    #    In addition to AOS and LOS subsatellite points, we calculate AOSbis and LOSbis, subsat points
    #    after certain detla seconds. This is used to calculate azimuth
    #
    # 2. calculate directions that are perpendicular (+90, -90 degrees) AZ_L, AZ_R
    #    (basic math, add/subtract 90 degrees, modulo 360)
    #
    # 3. calculate sensor swath (left-right "width" of the observation), divite by 2 to get D
    #    - SMAD
    #    - WERTZ Mission Geometry, page 420
    #
    # 4. calculate terminal distance starting from SSP at the azimuth AZ_L and AZ_R and distance D
    #    https://www.fcc.gov/media/radio/find-terminal-coordinates
    #    https://stackoverflow.com/questions/877524/calculating-coordinates-given-a-bearing-and-a-distance

    # TODO: Calculcate if this pass is northbound or southbound

    # Let's assume this is AVHRR instrument. Let's use its field of view angle.
    fov = AVHRR_FOV

    # Now calculate corner positions (use only the first method)
    swath = calc_swath(aos_lla[2], fov)
    print("Instrument angle is %f deg, altitude is %f km, swath (each side) is %f km, total swath is %f km" % (fov, aos_lla[2], swath, swath * 2))
    corner_ul = radial_distance(aos_lla[0], aos_lla[1], azimuth_add(aos_az, +90), swath)
    corner_ur = radial_distance(aos_lla[0], aos_lla[1], azimuth_add(aos_az, -90), swath)

    print("Upper left corner:  lat=%f lon=%f" % (corner_ul[0], corner_ul[1]))
    print("Upper right corner: lat=%f lon=%f" % (corner_ur[0], corner_ur[1]))

    # Now calculate corner positions (use only the first method)
    corner_ll = radial_distance(los_lla[0], los_lla[1], azimuth_add(los_az, +90), swath)
    corner_lr = radial_distance(los_lla[0], los_lla[1], azimuth_add(los_az, -90), swath)
    print("Lower left corner:  lat=%f lon=%f" % (corner_ll[0], corner_ll[1]))
    print("Lower right corner: lat=%f lon=%f" % (corner_lr[0], corner_lr[1]))

    # Ok, we have the sat position in LLA format. Getting sub-satellite point is trivial. Just assume altitude is 0.
    aos_lla = get_ssp(aos_lla)
    los_lla = get_ssp(los_lla)

    return d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr
