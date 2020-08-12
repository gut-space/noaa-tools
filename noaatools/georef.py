# This script processes images received from NOAA satellites

import sys
from datetime import datetime, timezone, timedelta
from math import atan, atan2, sqrt, pi, sin, cos, asin, acos, tan
from typing import Tuple

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.api import jday, Satrec

import cv2
from noaatools import imageproc

import numpy as np
from pymap3d import ecef

sys.path.append('.')
from noaatools import export_js
from noaatools.constants import DEG2RAD, RAD2DEG, Ellipsoid, Method, NOAA_PROCESSING_DELAY, RE, AVHRR_FOV, ellipsoid_wgs84

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
    T0 = 2451545.0 # J2000, 2000-Jan-01 12h UT1 as Julian date

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - T0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    T = d / 36525.0

    # Calculate GMST (in seconds at UT1=0)
    gmst = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 0.0000062 * T*T*T

    # Let's truncate this and return the value in degrees.
    # This is clearly broken.
    return (gmst % 24)*(15/3600.0)

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

    T0 = 2451545.0 # J2000, 2000-Jan-01 12h UT1 as Julian date

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - T0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    t = d / 36525.0

    # Don't understand this part.
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = ((jd + fr) % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau

    if theta > 2*pi:
        theta = theta - 2*pi
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
    return (lon + pi) % (2*pi) - pi

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

    lat = atan2(z, sqrt(x*x + y*y)) # phi
    lon = atan2(y, x) - gmst # lambda-E
    lon = longitude_trunc(lon)
    alt = sqrt(x*x + y*y + z*z) - RE # h

    # TODO: convert this to radians and use radians everywhere.
    return lat*RAD2DEG, lon*RAD2DEG, alt

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
    a  = ellipsoid.a
    f  = 1.0/ellipsoid.finv
    b  = a*(1 - 1.0/f)
    e2 = f*(2-f)

    phii = 1 # This is the starting value for initial iteration

    # There should be a check on |phii - phi| value, but let's always do 5 iterations. Good enough for now.
    for _ in range(1,5):

        C = 1/(sqrt(1-e2*pow(sin(phii), 2)))
        # This is not explicitly stated on clestrak page, but it's shown on a diagram.
        R = sqrt(x*x + y*y)
        phi = atan2(z + a*C*e2*sin(phii), R)
        h= R/(cos(phi)) - a*C

        phii=phi

    if type(t) == datetime:
        jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
        gmst = julianDateToGMST2(jd, fr)[0]
    elif type(t) == float:
        gmst = t
    else:
        raise ValueError("Incorrect parameter type of t")

    lon = atan2(y, x) - gmst # lambda-E
    lon = longitude_trunc(lon)

    return phi*RAD2DEG, lon*RAD2DEG, h

def teme2geodetic_pymap3d(x: float, y: float, z: float, t : datetime, ell = None):
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
    #lat, lon, alt = ecef.eci2geodetic(x*1000, y*1000, z*1000, t)
    #return lat, lon, alt/1000.0

    #print("teme[x,y,z]=%f, %f, %f" % (x, y, z))
    xecef, yecef, zecef = ecef.eci2ecef(np.array([x*1000]), np.array([y*1000]), np.array([z*1000]), t)
    #print("ecef[x,y,z]=%f, %f, %f" % (xecef, yecef, zecef))

    # True = we want the response in degrees
    lat, lon, alt = ecef.ecef2geodetic(xecef, yecef, zecef, ell, True)
    #print("lla = %f, %f, %f" % (lat, lon, alt))
    return lat, lon, alt/1000.0


def get_ssp(lla):
    return [ lla[0], lla[1], 0 ]

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

    # print("### calc_azimuth(lat1=%f lon1=%f, lat2=%f, lon2=%f)" % (lat1, lon1, lat2, lon2))

    d = 2*asin(sqrt((sin((lat1-lat2)/2))**2 +  cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))**2))

    tc1 = acos((sin(lat2)-sin(lat1)*cos(d))/(sin(d)*cos(lat1)))

    tc1 = atan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2))
    #if (tc1 < 0):
    #    tc1 += 2*pi
    if (tc1 > 2*pi):
        tc1 -= 2*pi
    return tc1*RAD2DEG

def calc_swath(alt, nu):
    """This calculates the swath width, given the altitude (alt, in km) of the sat and camera angle (nu, in radians).
        Returns swath in km"""

    # Convert to radians first.
    nur = nu*DEG2RAD

    # Ok, this is an overly simplified approximation. It neglects the Earth curvature.
    # return alt*tan(nu)

    # Source Wertz "Mission geometry", pg. 420.

    # rho is an angle between two lines: (sat - tangential to Earth) and (sat - Earth center)
    rho = asin(RE/(RE+alt))

    epsilon = acos(sin(nur)/sin(rho))
    lam = pi/2 - nur - epsilon
    swath = RE*lam
    print("calc_swath(alt=%f nu=%f/%f) => rho=%f/%f epsilon=%f/%f, lambda= %f/%f => swath=%f [km]" %
    (alt, nu, nur, rho, rho*RAD2DEG, epsilon, epsilon*RAD2DEG, lam, lam*RAD2DEG, swath))

    return swath

def radial_distance(lat1, lon1, bearing, distance):
    """
    Return final coordinates (lat2,lon2) [in degrees] given initial coordinates
    (lat1,lon1) [in degrees] and a bearing [in degrees] and distance [in km]

    Based on this:
    https://stackoverflow.com/questions/877524/calculating-coordinates-given-a-bearing-and-a-distance
    """

    rlat1 = lat1*DEG2RAD
    rlon1 = lon1*DEG2RAD
    rdistance = distance / RE # normalize linear distance to radian angle
    rbearing = bearing * DEG2RAD

    rlat = asin( sin(rlat1) * cos(rdistance) + cos(rlat1) * sin(rdistance) * cos(rbearing) )

    if cos(rlat) == 0 or abs(cos(rlat)) < 0.00000001: # Endpoint a pole
        rlon = rlon1
    else:
        rlon = ( (rlon1 + asin( sin(rbearing)* sin(rdistance) / cos(rlat) ) + pi ) % (2*pi) ) - pi

    return (rlat*RAD2DEG, rlon*RAD2DEG)

def calc_distance(lat1, lon1, lat2, lon2):
    """
    Calculates distance between two (lat,lon) points, expressed in rad. Return value is in km.
    """
    rlat1 = lat1*DEG2RAD
    rlon1 = lon1*DEG2RAD
    rlat2 = lat2*DEG2RAD
    rlon2 = lon2*DEG2RAD

    d = 2 * asin(sqrt((sin((rlat1-rlat2)/2))**2 + cos(rlat1)*cos(rlat2)*(sin((rlon1-rlon2)/2))**2))

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

def azimuth_apt(lat1, lon1, lat2, lon2):
    """everythin in rad """
    delta_lon = lon2 - lon1

    az = sin(delta_lon) / ( cos(lat1)*tan(lat2) - sin(lat1)*cos(delta_lon) )

    return az

def distance_apt(lat1, lon1, lat2, lon2):

    delta_lon = lon2 - lon1

    cos_central_angle = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(delta_lon)

    cos_central_angle = max(cos_central_angle, -1)
    cos_central_angle = min(cos_central_angle, 1)

    return acos(cos_central_angle)


def latlon_to_rel_px(latlon, start_latlon, ref_az, xres, yres, yaw) -> (float, float):
    """
    latlon - a tuple of lat, lon, alt
    start_latlon - a tuple of lat, lon, alt
    ref_az - reference azimuth (in radians)
    xres - horizontal resolution per pixel
    yres - vertical resolution per pixel
    """

    az = calc_azimuth(start_latlon, latlon) * DEG2RAD
    B = az - ref_az

    # TODO: there's this weird shit in noaa-apt: map.rs:106
    #
    #    // Set maximum, otherwise we get wrapping problems I do not fully
    #    // understand: opposite parts of the world are mapped to the same
    #    // position because of the cyclic nature of sin(), cos(), etc.
    #    let c = geo::distance(latlon, start_latlon).max(-PI/3.).min(PI/3.);
    c1 = calc_distance(latlon[0]*RAD2DEG, latlon[1]*RAD2DEG, start_latlon[0]*RAD2DEG, start_latlon[1]*RAD2DEG) / RE
    c11= distance_apt(latlon[0]*RAD2DEG, latlon[1]*RAD2DEG, start_latlon[0]*RAD2DEG, start_latlon[1]*RAD2DEG)
    c2 = max(c1, -pi/3.0)
    c = min(c2, pi/3.0)
    print("### distance(latlon=(%f, %f), start_latlon=(%f, %f)) c1=%f c11=%f c2=%f c=%f" % (latlon[0], latlon[1], start_latlon[0], start_latlon[1], c1, c11, c2, c))

    a = atan(cos(B) * tan(c))
    b = asin(sin(B) * sin(c))

    x = -b / xres
    y = a / yres + yaw * x

    print("### latlon_to_rel_px: ref_az=%f, az=%f, B=%f c=%f a=%f b=%f x=%f y=%f" % (ref_az, az, B, c, a, b, x, y))

    return x,y

def draw_line(image, latlon1, latlon2, rgba, ref_az, xres, yres, yaw, sat_positions):

    start_latlon = sat_positions[0]
    print("##### start_latlon=%f,%f" % (start_latlon[0], start_latlon[1]))

    # Convert latlon to (x, y) to pixel coordinates
    x1, y1 = latlon_to_rel_px(latlon1, start_latlon, ref_az, xres, yres, yaw)
    x2, y2 = latlon_to_rel_px(latlon2, start_latlon, ref_az, xres, yres, yaw)

    print("## draw_line: latlon1=(%f,%f) latlon2=(%f,%f), x1=%f y1=%f, x2=%f, y2=%f" % (latlon1[0], latlon1[1], latlon2[0], latlon2[1], x1, y1, x2, y2))

    height, width, _ = image.shape

    # Offset correction on y
    est_y1 = int (min( max(y1, 0.), height - 1)) # make sue y1 is in <0...height-1> range
    est_y2 = int (min( max(y2, 0.), height - 1))

    print("### draw_line() est_y1=%d, est_y2=%d" % (est_y1, est_y2))

    x1_offset, _ = latlon_to_rel_px(sat_positions[est_y1], start_latlon, ref_az, xres, yres, yaw)
    x2_offset, _ = latlon_to_rel_px(sat_positions[est_y2], start_latlon, ref_az, xres, yres, yaw)
    x1 -= x1_offset
    x2 -= x2_offset

    # See if at least one point is inside
    if (x1 > -456 and x1 < 456 and y1 > 0. and y1 < height) or (x1 > -600. and x1 < 600. and y1 > 0. and y1 < height):

        cv2.line(image, (x1 + 539,y1), (x2 + 539,y2) , (0,0,255), 5)
        cv2.line(image, (x1 + 1579,y1), (x2 + 1579,y2) , (0,0,255), 5)


def georef_apt(method: Method, tle1: str, tle2: str, aos_txt: str, los_txt: str, imgfile: str):
    """ This georeferencing method is roughly based on noaa-apt, see
        https://noaa-apt.mbernardi.com.ar/"""

    img = cv2.imread(imgfile)
    height, width, _ = img.shape

    # Convert date as a string datetime. Make sure to use UTC rather than the default (local timezone)
    d1 = datetime.fromisoformat(aos_txt).replace(tzinfo=timezone.utc)

    # This approach uses new API 2.x which gives a slightly different results.
    # In case of NOAA, the position is off by less than milimeter
    sat = Satrec.twoline2rv(tle1, tle2)
    jd, fr = jday(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)

    print("## AOS=%s" % d1)

    sat_positions = []

    print("## There are %d lines in the %s file." % (height, imgfile))
    for i in range(0,height):
        _, pos, _ = sat.sgp4(jd, fr)
        gmst, _ = julianDateToGMST2(jd, fr)
        lla = teme2geodetic_oblate(pos[0], pos[1], pos[2], gmst, ellipsoid_wgs84)
        lla = (lla[0]*DEG2RAD, lla[1]*DEG2RAD, lla[2])
        if i % 10 == 11:
            print("## %d of %d: t=%s result=%s gmst=%f lon=%f lat=%f alt=%f" % (i, height, d1, pos, gmst, lla[0], lla[1], lla[2]))
        # TODO: print timestamp
        sat_positions.append(lla) # degrees
        fr = fr + 0.5/86400.0

    rows = len(sat_positions)
    print("## %d positions calculated" % rows)
    print("## Sat position for row 0: %f, %f" % (sat_positions[0][0], sat_positions[0][1])) # degrees

    print("## Sat position for row %d: %f, %f" % (rows, sat_positions[rows -1 ][0], sat_positions[rows - 1][1]))

    ref_az = calc_azimuth(sat_positions[0], sat_positions[-1])*DEG2RAD
    if (ref_az < 0):
        ref_az += 2*pi

    dist = calc_distance(sat_positions[0][0], sat_positions[0][1], sat_positions[-1][0], sat_positions[-1][1])
    dist /= RE

    print("## azimuth=%f, dist=%f height=%f" % (ref_az / 180*pi, dist, height))

    vscale = 1.0000 # vertical scale
    hscale = 1.0000 # horizontal scale

    yres = dist / float(height) / vscale
    xres = 0.0005 / hscale
    yaw = 0.0

    print("## xres=%f, yres=%f, yaw=%f" % (xres, yres, yaw))

    print("------------")

    img_mod = draw_line(img, (53.91*DEG2RAD, 14.28*DEG2RAD), (44.92*DEG2RAD, 12.55*DEG2RAD), (0,0,255),
                        ref_az, xres, yres, yaw, sat_positions)

    cv2.imshow("Line",img_mod)
    cv2.waitKey(0)
    cv2.destroyAllWindows()


def georef(method: Method, tle1: str, tle2: str, aos_txt: str, los_txt: str, imgfile: str):
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
    #sat_old = twoline2rv(tle1, tle2, wgs72)
    #pos1_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    #pos2_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)

    # This approach uses new API 2.x which gives a slightly different results.
    # In case of NOAA, the position is off by less than milimeter
    sat = Satrec.twoline2rv(tle1, tle2)
    jd1, fr1 = jday(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    jd2, fr2 = jday(d2.year, d2.month, d2.day, d2.hour, d2.minute, d2.second)

    # Take sat processing/transmission delay into consideration. At AOS time the signal received
    # was already NOAA_PROCESSING_DELAY seconds old.
    fr1 -= NOAA_PROCESSING_DELAY/86400.0
    fr2 -= NOAA_PROCESSING_DELAY/86400.0

    _, pos1, _ = sat.sgp4(jd1, fr1) # returns error, position and velocity - we care about position only
    _, pos2, _ = sat.sgp4(jd2, fr2)

    # Delta between a point and a point+delta (the second delta point is used to calculate azimuth)
    DELTA = 30.0

    _, pos1delta, _ = sat.sgp4(jd1, fr1 + DELTA/86400.0)
    _, pos2delta, _ = sat.sgp4(jd2, fr2 + DELTA/86400.0)

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

    # Ok, we have sat position at time of AOS and LOS returned by SGP4 models. The tricky part here is those are in
    # Earth-Centered Intertial (ECI) reference system. The model used is TEME (True equator mean equinox).

    # Convert AOS coordinates to LLA
    aos_lla = teme2geodetic(method, pos1[0], pos1[1], pos1[2], d1)

    # Now caluclate a position for AOS + 30s. Will use it to determine the azimuth
    d1delta = d1 + timedelta(seconds = 30.0)
    aos_bis = teme2geodetic(method, pos1delta[0], pos1delta[1], pos1delta[2], d1delta)
    aos_az = calc_azimuth(aos_lla, aos_bis)

    print("AOS converted to LLA is lat=%f long=%f alt=%f, azimuth=%f" % (aos_lla[0], aos_lla[1], aos_lla[2], aos_az) )

    # Now do the same for LOS
    los_lla = teme2geodetic(method, pos2[0], pos2[1], pos2[2], d2)

    # Let's use a point 30 seconds later. AOS and (AOS + 30s) will determine the azimuth
    d2delta = d2 + timedelta(seconds = 30.0)
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
    print("Instrument angle is %f deg, altitude is %f km, swath (each side) is %f km, total swath is %f km" % (fov, aos_lla[2], swath, swath*2))
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
