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
from datetime import datetime, timezone
from collections import namedtuple
from math import atan, atan2, sqrt, pi, sin, cos

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.api import jday, Satrec

# This defines ellipsoid (a = equatorial radius in km, finv = inverse flattening)
Ellipsoid = namedtuple('Ellipsoid', "a finv")

# Source: https://en.wikipedia.org/wiki/Earth_ellipsoid#Historical_Earth_ellipsoids
ellipsoid_wgs84 = Ellipsoid(a = 6378.137, finv = 298.257223563)

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

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - 2451545.0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    T = d / 36525

    # Calculate GMST (in seconds at UT1=0)
    gmst = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 0.0000062 * T*T*T

    # Let's truncate this and return the value in degrees.
    return (gmst % 24)*(15/3600.0)

def temeToGeodetic_spherical(x, y, z, jd, fr):
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

    gmst = julianDateToGMST(jd, fr)

    RE = 6378.137 # Earth radius (in km)

    lat = atan2(z, sqrt(x*x + y*y)) # phi
    lon = atan2(y, x) - gmst # lambda-E
    alt = sqrt(x*x + y*y + z*z) - RE # h

    return lat*180/pi, lon*180/pi, alt

def temeToGeodetic(x, y, z, ellipsoid, jd, fr):
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

    gmst = julianDateToGMST(jd, fr)

    lon = atan2(y, x) - gmst # lambda-E

    return phi*180/pi, lon*180/pi, h

def georef_naive(tle1, tle2, aos, los):
    """ This is a naive georeferencing method:
        - calculates the sat location at AOS and LOS points (using )
    then calculates distance between them. """

    d1 = datetime.fromisoformat(aos).replace(tzinfo=timezone.utc)
    d2 = datetime.fromisoformat(los).replace(tzinfo=timezone.utc)

    # Sat image aquisition and image transmission is not instanteneous. There is some delay. It's small, but
    # bacause the sat moves at speeds of more than 7km/s, even fractions of second make a big deal of a difference.
    # Not sure if this data is available anywhere. I think it'll have to be picked with trial-and-error.
    # This is expressed in seconds.
    delay = 0.5

    # STEP 1: Calculate sat location at AOS and LOS

    # This approach uses old API 1.x
    #sat_old = twoline2rv(tle1, tle2, wgs72)
    #pos1_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    #pos2_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)

    # This approach uses new API 2.x which gives a slightly different results.
    # In case of NOAA, the position is off by less than milimeter
    sat = Satrec.twoline2rv(tle1, tle2)
    jd1, fr1 = jday(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    jd2, fr2 = jday(d2.year, d2.month, d2.day, d2.hour, d2.minute, d2.second)
    fr1 += delay/86400.0 # Take sat processing/transmission delay into consideration
    fr2 += delay/86400.0

    _, pos1, _ = sat.sgp4(jd1, fr1) # returns error, position and velocity - we care about position only
    _, pos2, _ = sat.sgp4(jd2, fr2)

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

    # Ok, we have sat position at time of AOS and LOS. The tricky part here is those are in
    # Earth-Centered Intertial (ECI) reference system. The model used is TEME (True equator,
    # mean equinox). Let's calculate the SSP using two methods:
    # - simple (that assumes spherical Earth)
    # - oblate Earth (uses passed ellipsoid, WGS84 in this case)

    print("METHOD 1 (spherical Earth)")
    lla1 = temeToGeodetic_spherical(pos1[0], pos1[1], pos1[2], jd1, fr1)
    print("AOS: ECI[x=%f, y=%f, z=%f] converted to LLA is long=%f lat=%f alt=%f" %
    (pos1[0], pos1[1], pos1[2], lla1[0], lla1[1], lla1[2]))

    lla2 = temeToGeodetic_spherical(pos2[0], pos2[1], pos2[2], jd2, fr2)
    print("LOS: ECI[x=%f, y=%f, z=%f] converted to LLA is long=%f lat=%f alt=%f" %
    (pos2[0], pos2[1], pos2[2], lla2[0], lla2[1], lla2[2]))

    print("METHOD 2 (oblate Earth)")
    lla1 = temeToGeodetic(pos1[0], pos1[1], pos1[2], ellipsoid_wgs84, jd1, fr1)
    print("AOS: ECI[x=%f, y=%f, z=%f] converted to LLA is long=%f lat=%f alt=%f" %
    (pos1[0], pos1[1], pos1[2], lla1[0], lla1[1], lla1[2]))

    lla2 = temeToGeodetic(pos2[0], pos2[1], pos2[2], ellipsoid_wgs84, jd2, fr2)
    print("LOS: ECI[x=%f, y=%f, z=%f] converted to LLA is long=%f lat=%f alt=%f" %
    (pos2[0], pos2[1], pos2[2], lla2[0], lla2[1], lla2[2]))

    # STEP 3: Calculate the radial distance between AOS SSP and LOS SSP, divide is by image height. The result will be
    # angular resolution per pixel. Now multiply the value by image width/2 and then add/subtract from the AOS/LOS SSP
    # to get corners of the image.

    # STEP 4: Export georeferencing data.

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

    if len(sys.argv) < 5:
        usage()
        sys.exit(-1)

    #process(sys.argv[1], params)

    tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
    tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

    aos = '2020-04-12 09:01:03.063476'
    los = '2020-04-12 09:17:06.466954'

    georef_naive(tle1, tle2, aos, los)


# NOTES: If the above doesn't work, the next thing to try will be this:
# https://github.com/geospace-code/pymap3d, especially eci2geodetic

# This one is interesting.
# https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_geodetic_coordinates


# This one has almost complete solution using Skyfield
# https://space.stackexchange.com/questions/19339/better-way-to-get-approximate-ground-track-for-a-satellite-using-skyfield/40218#40218
