# This script processes images received from NOAA satellites
USAGE = '''
#
# noaa_process script
#
# Current capabilities:
# - basic sanity checking (if the image is indeed NOAA observation)
# - mark left and right images
# - extract left and right images
# - perform histogram stretch
# - display images
#
# Usage:
# python noaa_process.py filename.png tle.txt AOS LOS
#
# filename.png - this is a file that's decoded with noaa-apt
# tle.txt - TLE information of the orbit
# aos - aquisition of signal
# los - loss of signal
'''

import sys
from datetime import datetime, timezone, timedelta
from collections import namedtuple
from math import atan, atan2, sqrt, pi, sin, cos, asin, acos
from enum import Enum

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.api import jday, Satrec

import numpy as np
from pymap3d import ecef

from export_czml import *

class Methods(Enum):
    SPHERICAL = 1
    OBLATE = 2
    PYMAP3D = 3

# Conversion between radians and degrees
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321

RE = 6378.137 # Earth radius (in km)

# This defines ellipsoid (a = equatorial radius in km, finv = inverse flattening)
Ellipsoid = namedtuple('Ellipsoid', "a finv")

# Source: https://en.wikipedia.org/wiki/Earth_ellipsoid#Historical_Earth_ellipsoids
ellipsoid_wgs84 = Ellipsoid(a = RE, finv = 298.257223563)

# Sat image aquisition and image transmission is not instanteneous. There is some delay. It's small, but
# bacause the sat moves at speeds of more than 7km/s, even fractions of second make a big deal of a difference.
# Not sure if this data is available anywhere. I think it'll have to be picked with trial-and-error.
# This is expressed in seconds.
NOAA_PROCESSING_DELAY = 0.0

#######################################################################################################
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

def julianDateToGMST2(jd, fr):
    """
    Converts Julian date (expressed at two floats) to GMST (Greenwich Mean Sidereal Time 1982).

    Parameters:
    jd : float - Julian date full integer + 0.5
    fr : float - fractional part of the Julian date

    Returns
    =======
    A single floating point representing a GMST, expressed in radians

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

    # Don't undersran
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t
    dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t
    theta = ((jd + fr) % 1.0 + g * _second % 1.0) * tau
    theta_dot = (1.0 + dg * _second / 36525.0) * tau
    return theta, theta_dot

def longitude_trunc(lon):
    """ Makes sure the longitude is within -2*pi ... 2*pi range """
    if (lon < -pi):
        lon += 2*pi
    if (lon > pi):
        lon -= 2*pi
    return lon

def teme2geodetic_spherical(x, y, z, t):
    """
    Converts ECI coords (x,y,z - expressed in km) to LLA (longitude, lattitude, altitude).
    This function assumes the Earth is completely round.

    The calculations here are based on T.S. Kelso's excellent paper "Orbital Coordinate Systems, Part III
    https://celestrak.com/columns/v02n03/.

    Parameters
    ==========
    x,y,z : float - coordates in TEME (True Equator Mean Equinoex) version of ECI (Earth Centered Intertial) coords system.
            This is the system that's produced by SGP4 models.
    jd, fr : float - julian date - expressed as two floats that should be summed together.
    """

    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    gmst = julianDateToGMST2(jd, fr)[0]

    lat = atan2(z, sqrt(x*x + y*y)) # phi
    lon = atan2(y, x) - gmst # lambda-E
    lon = longitude_trunc(lon)
    alt = sqrt(x*x + y*y + z*z) - RE # h

    return lat*RAD2DEG, lon*RAD2DEG, alt

def teme2geodetic_oblate(x, y, z, t, ellipsoid):
    """
    Converts ECEF coords (x,y,z - expressed in km) to LLA (longitude, lattitude, altitude).
    ellipsoid is Earth ellipsoid to be used (e.g. ellipsoid_wgs84).

    The calculations here are based on T.S. Kelso's excellent paper "Orbital Coordinate Systems, Part III
    https://celestrak.com/columns/v02n03/.

    Parameters
    ==========
    x,y,z : float - coordates in TEME (True Equator Mean Equinoex) version of ECI (Earth Centered Intertial) coords system.
            This is the system that's produced by SGP4 models.
    ellipsoid: Ellipsoid - an Earth exlipsoid specifying Earth oblateness, e.g. Ellipsoid_wgs84. Two params are used from it:
            a and inverse of f. Both must be specified in kms
    jd, fr : float - julian date - expressed as two floats that should be summed together.

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

    # There should be a check on |phii - phi| value, but let's always do 4 iterations. Good enough for now.
    for iter in range(1,5):

        C = 1/(sqrt(1-e2*pow(sin(phii), 2)))
        # This is not explicitly stated on clestrak page, but it's shown on a diagram.
        R = sqrt(x*x + y*y)
        phi = atan2(z + a*C*e2*sin(phii), R)
        h= R/(cos(phi)) - a*C

        # print("iter=%d, phi=%f, phii=%f, |phi-phii|=%f" % (iter, phi, phii, abs(phi-phii)))

        phii=phi

    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    gmst = julianDateToGMST2(jd, fr)[0]

    lon = atan2(y, x) - gmst # lambda-E
    lon = longitude_trunc(lon)

    return phi*RAD2DEG, lon*RAD2DEG, h

def teme2geodetic_pymap3d(x, y, z, t : datetime, ell = None):
    """
    Converts TEME coordinates to geodetic, using pymap3d library.
    For details, see https://github.com/geospace-code/pymap3d

    Returns
    =======
    lat, lon, alt - latitude, longitude (both in degrees), alt (in km)
    """

    # Short version - whole conversion in one go
    #lat, lon, alt = ecef.eci2geodetic(x*1000, y*1000, z*1000, t)
    #return lat, lon, alt/1000.0

    #print("teme[x,y,z]=%f, %f, %f" % (x, y, z))
    xecef, yecef, zecef = ecef.eci2ecef(np.array([x*1000]), np.array([y*1000]), np.array([z*1000]), t)
    #print("exef[x,y,z]=%f, %f, %f" % (xecef, yecef, zecef))

    # True = we want the response in degrees
    lat, lon, alt = ecef.ecef2geodetic(xecef, yecef, zecef, ell, True)
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

    d = 2*asin(sqrt((sin((lat1-lat2)/2))**2 +  cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))**2))

    tc1 = acos((sin(lat2)-sin(lat1)*cos(d))/(sin(d)*cos(lat1)))

    tc1 = atan2(sin(lon1-lon2)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2))
    if (tc1 < 0):
        tc1 += 2*pi
    if (tc1 > 2*pi):
        tc1 -= 2*pi
    return tc1*RAD2DEG

def calculate_swath(alt, deg):
    """This calculates the swath width, given the altitude of the sat and camera angle.
        Returns swath in km"""

    # TODO
    return 200

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

def azimuth_add(az, delta):
    """ Adds delta to specified azimuth. Does the modulo 360 arithmetic"""

    return (az + delta) % 360.0

def georef(imgname, tle1, tle2, aos, los):
    """ This is a naive georeferencing method:
        - calculates the sat location at AOS and LOS points (using )
    then calculates distance between them. """

    # Convert date as a string datetime. Make sure to use UTC rather than the default (local timezone)
    d1 = datetime.fromisoformat(aos).replace(tzinfo=timezone.utc)
    d2 = datetime.fromisoformat(los).replace(tzinfo=timezone.utc)

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

    _, pos1delta, _ = sat.sgp4(jd1, fr1 + DELTA/86400)
    _, pos2delta, _ = sat.sgp4(jd2, fr2 + DELTA/86400)

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

    # Ok, we have sat position at time of AOS and LOS. The tricky part here is those are in
    # Earth-Centered Intertial (ECI) reference system. The model used is TEME (True equator,
    # mean equinox). Let's calculate the SSP using three methods:
    # - simple (that assumes spherical Earth)
    # - oblate Earth (uses passed ellipsoid, WGS84 in this case)
    # - using pymap3d lib (which is most precise)

    # AOS calc
    print("AOS: ECI[x=%f, y=%f, z=%f]" % (pos1[0], pos1[1], pos1[2]))

    aos1 = teme2geodetic_spherical(pos1[0], pos1[1], pos1[2], d1)
    aos2 = teme2geodetic_oblate(pos1[0], pos1[1], pos1[2], d1, ellipsoid_wgs84)
    aos3 = teme2geodetic_pymap3d(pos1[0], pos1[1], pos1[2], d1)

    print("METHOD 1 (spherical Earth): AOS converted to LLA is lat=%f long=%f alt=%f" % (aos1[0], aos1[1], aos1[2]) )
    print("METHOD 2 (oblate Earth):    AOS converted to LLA is lat=%f long=%f alt=%f" % (aos2[0], aos2[1], aos2[2]) )
    print("METHOD 3 (pymap3d):         AOS converted to LLA is lat=%f long=%f alt=%f" % (aos3[0], aos3[1], aos3[2]) )

    # Let's use a point 30 seconds later. AOS and (AOS + 30s) will determine the azimuth
    d1delta = d1 + timedelta(seconds = 30.0)

    # Now calculate azimuth
    aos1bis = teme2geodetic_spherical(pos1delta[0], pos1delta[1], pos1delta[2], d1delta)
    aos2bis = teme2geodetic_oblate(pos1delta[0], pos1delta[1], pos1delta[2], d1delta, ellipsoid_wgs84)
    aos3bis = teme2geodetic_pymap3d(pos1delta[0], pos1delta[1], pos1delta[2], d1delta)
    az1 = calc_azimuth(aos1, aos1bis)
    az2 = calc_azimuth(aos2, aos2bis)
    az3 = calc_azimuth(aos3, aos3bis)

    print("METHOD1 AOS azimuth = %f" % az1)
    print("METHOD2 AOS azimuth = %f" % az2)
    print("METHOD3 AOS azimuth = %f" % az3)

    avhrr_angle = 54.3

    # Now calculate corner positions (use only the first method)
    corner_ul = radial_distance(aos1[0], aos1[1], az1 + 90, calculate_swath(aos1[2], avhrr_angle / 2.0 ))
    corner_ur = radial_distance(aos1[0], aos1[1], az1 - 90, calculate_swath(aos1[2], avhrr_angle / 2.0 ))

    print("METHOD1 Upper left corner:  lat=%f lon=%f" % (corner_ul[0], corner_ul[1]))
    print("METHOD1 Upper right corner: lat=%f lon=%f" % (corner_ur[0], corner_ur[1]))

    # LOS
    los1 = teme2geodetic_spherical(pos2[0], pos2[1], pos2[2], d2)
    los2 = teme2geodetic_oblate(pos2[0], pos2[1], pos2[2], d2, ellipsoid_wgs84)
    los3 = teme2geodetic_pymap3d(pos2[0], pos2[1], pos2[2], d2)

    print("METHOD 1 (spherical Earth): LOS converted to LLA is lat=%f long=%f alt=%f" % (los1[0], los1[1], los1[2]))
    print("METHOD 2 (oblate Earth):    LOS converted to LLA is lat=%f long=%f alt=%f" % (los2[0], los2[1], los2[2]))
    print("METHOD 3 (pymap3d):         LOS converted to LLA is lat=%f long=%f alt=%f" % (los3[0], los3[1], los3[2]))

    # Let's use a point 30 seconds later. AOS and (AOS + 30s) will determine the azimuth
    d2delta = d2 + timedelta(seconds = 30.0)

    los1bis = teme2geodetic_spherical(pos2delta[0], pos2delta[1], pos2delta[2], d2delta)
    los2bis = teme2geodetic_oblate(pos2delta[0], pos2delta[1], pos2delta[2], d2delta, ellipsoid_wgs84)
    los3bis = teme2geodetic_pymap3d(pos2delta[0], pos2delta[1], pos2delta[2], d2delta)
    print("METHOD 1 (spherical Earth): LOSbis converted to LLA is lat=%f long=%f alt=%f" % (los1bis[0], los1bis[1], los1bis[2]))
    print("METHOD 2 (oblate Earth):    LOSbis converted to LLA is lat=%f long=%f alt=%f" % (los2bis[0], los2bis[1], los2bis[2]))
    print("METHOD 3 (pymap3d):         LOSbis converted to LLA is lat=%f long=%f alt=%f" % (los3bis[0], los3bis[1], los3bis[2]))

    az1los = calc_azimuth(los1, los1bis)
    az2los = calc_azimuth(los2, los2bis)
    az3los = calc_azimuth(los3, los3bis)

    print("METHOD1 LOS azimuth = %f" % az1los)
    print("METHOD2 LOS azimuth = %f" % az2los)
    print("METHOD3 LOS azimuth = %f" % az3los)

    # Now calculate corner positions (use only the first method)
    corner_ll = radial_distance(los1[0], los1[1], az1los + 90, calculate_swath(los1[2], avhrr_angle / 2.0 ))
    corner_lr = radial_distance(los1[0], los1[1], az1los - 90, calculate_swath(los1[2], avhrr_angle / 2.0 ))
    print("METHOD1 Lower left corner:  lat=%f lon=%f" % (corner_ll[0], corner_ll[1]))
    print("METHOD1 Lower right corner: lat=%f lon=%f" % (corner_lr[0], corner_lr[1]))

    # Ok, we have the sat position in LLA format. Getting sub-satellite point is trivial. Just assume altitude is 0.
    aos1 = get_ssp(aos1)
    aos2 = get_ssp(aos2)
    aos3 = get_ssp(aos3)

    los1 = get_ssp(los1)
    los2 = get_ssp(los2)
    los3 = get_ssp(los3)

    # STEP 3: Find image corners. Here's an algorithm proposal:
    #
    # 1. Calculate the orthodromic distance (a great-circle distance) between AOS SSP and LOS SSP.
    #    See https://en.wikipedia.org/wiki/Great-circle_distance
    #
    # 2. Find great circle that passes through AOS SSP and LOS SSP.
    #
    # 3. Find a great circle that passes through AOS SSP and is perpendicular to great circle found in step 2.
    #    Calculate cross-track distance (use image's widht/height multiplied by the orthodromic distance in step 1),
    #    traverse half the distance in each way to find two corners of the image.
    #
    # 4. Repreat step 3, but use LOS SSP. This will find the remaining two corners of the image.
    #
    # Alternative algorithm:
    # 1. calculate satellite flight azimuth AZ
    #    https://en.wikipedia.org/wiki/Great-circle_navigation
    # 2. calculate directions that are perpendicular (+90, -90 degrees) AZ_L, AZ_R
    #    (basic math, add/subtract 90 degrees, modulo 360)
    # 3. calculate sensor swath (left-right "width" of the observation), divite by 2 to get D
    #    - SMAD
    #    - WERTZ Mission Geometry
    # 4. calculate terminal distance starting from SSP at the azimuth AZ_L and AZ_R and distance D
    #    https://www.fcc.gov/media/radio/find-terminal-coordinates
    #    https://stackoverflow.com/questions/877524/calculating-coordinates-given-a-bearing-and-a-distance

    # STEP 4: Do we need to detect if it's northbound or southbound fly-over? We can do it easily, by
    # comparing aos_lat (lla1[1]) with los_lat (lla2[1])

    # STEP 5: Export georeferencing data.
    outfile = ".".join(imgname.split('.')[:-1]) + ".js"
    aos_list = [ aos1, aos2, aos3 ]
    los_list = [ los1, los2, los3 ]
    methods  = [ "spherical", "oblate", "pymap3d" ]
    export2cesium(outfile, imgname, d1, d2, aos_list, los_list, methods, tle1, tle2, corner_ul, corner_ur, corner_ll, corner_lr)

    # STEP 6: (possibly outside of this script):
    # - use GDAL library to georeference image (https://pcjericks.github.io/py-gdalogr-cookbook/)
    # - import into QGIS (and follow this tutorial: https://www.qgistutorials.com/en/docs/georeferencing_basics.html)
    # - display the image in Cesium (https://www.cesium.com/docs/cesiumjs-ref-doc/SingleTileImageryProvider.html)

def usage():
    print(USAGE)

if __name__ == "__main__":
    params = {
        "histogram": True,
        "histogram-adaptive": True,
        "border": True,
        "show": True,
        "write": False,
        "write-left": True,
        "write-right": False,
        "denoise": False,
        "georef": True # Georeference
    }

    # Let's ignore input parameters and pretend we were asked to georeference observation #1276.
    #if len(sys.argv) < 5:
    #    usage()
    #    sys.exit(-1)

    tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
    tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

    aos = '2020-04-12 09:01:03.063476'
    los = '2020-04-12 09:17:06.466954'

    georef("1276.png", tle1, tle2, aos, los)


# This one is interesting.
# https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_geodetic_coordinates


# This one has almost complete solution using Skyfield
# https://space.stackexchange.com/questions/19339/better-way-to-get-approximate-ground-track-for-a-satellite-using-skyfield/40218#40218
