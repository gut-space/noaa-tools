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

import numpy as np
from pymap3d import ecef

# This is needed to export orbit to CZML format (Cesium)
from tletools import TLE
from poliastro.czml.extract_czml import CZMLExtractor
from astropy import time

# This defines ellipsoid (a = equatorial radius in km, finv = inverse flattening)
Ellipsoid = namedtuple('Ellipsoid', "a finv")

# Source: https://en.wikipedia.org/wiki/Earth_ellipsoid#Historical_Earth_ellipsoids
ellipsoid_wgs84 = Ellipsoid(a = 6378.137, finv = 298.257223563)

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

    # First calculate number of days since J2000 (2000-Jan-01 12h UT1)
    d = jd - 2451545.0
    d = d + fr

    # Now convert this to centuries. Don't ask me why.
    T = d / 36525

    # Calculate GMST (in seconds at UT1=0)
    gmst = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 0.0000062 * T*T*T

    # Let's truncate this and return the value in degrees.
    return (gmst % 24)*(15/3600.0)

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
    gmst = julianDateToGMST(jd, fr)

    RE = 6378.137 # Earth radius (in km)

    lat = atan2(z, sqrt(x*x + y*y)) # phi
    lon = atan2(y, x) - gmst # lambda-E
    alt = sqrt(x*x + y*y + z*z) - RE # h

    return lat*180/pi, lon*180/pi, alt

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
    gmst = julianDateToGMST(jd, fr)

    lon = atan2(y, x) - gmst # lambda-E

    return phi*180/pi, lon*180/pi, h

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

    print("teme[x,y,z]=%f, %f, %f" % (x, y, z))
    xecef, yecef, zecef = ecef.eci2ecef(np.array([x*1000]), np.array([y*1000]), np.array([z*1000]), t)
    print("exef[x,y,z]=%f, %f, %f" % (xecef, yecef, zecef))

    # True = we want the response in degrees
    lat, lon, alt = ecef.ecef2geodetic(xecef, yecef, zecef, ell, True)
    return lat, lon, alt/1000.0

def cesium_preamble():
    code = """
    // Code generated with noaa-tools.
    // This can be executed using Cesium. The easiest way to test it is to go to
    // https://sandcastle.cesium.com and paste this code there.

    var viewer = new Cesium.Viewer('cesiumContainer');
    var pinBuilder = new Cesium.PinBuilder();
    """

    return code

def export2cesium_point(lla, name):
    """
    Exports specified LLA coordinates to filename, using name as a label.
    """

    code = """
    var questionPin = viewer.entities.add({
        name : '%s',
        position : Cesium.Cartesian3.fromDegrees(%f, %f, %f),
        billboard : {
            image : pinBuilder.fromText('%s', Cesium.Color.RED, 48).toDataURL(),
            verticalOrigin : Cesium.VerticalOrigin.BOTTOM
        }
    });
    """ % (name, lla[0], lla[1], lla[2], name[0])

    return code

def export2cesium_tle(tle1, tle2, satname, aos, los):
    """
    Generates JS Cesium code that shows satellite trajectory between AOS and LOS.

    Parameters
    ==========
    tle1, tle2 - first and second line of TLE
    satname - to be displayed
    aos - Aquisition of Singal in Timedate format
    los - Loss of Signal in Timedate format
    """
    tle = TLE.from_lines(satname, tle1, tle2)
    orb = tle.to_orbit()

    sample_points = 10

    aos_astropy = time.Time(aos, scale="utc")
    los_astropy = time.Time(los, scale="utc")

    extractor = CZMLExtractor(aos_astropy, los_astropy, sample_points)
    extractor.add_orbit(orb, path_show=True, path_width=3, path_color=[125, 80, 120, 255], label_text=satname)

    txt =" var czml = [\n"
    for i in extractor.packets:
        txt += repr(i)
        txt += ",\n"
    txt += "];\n"

    # I don't understand what's exactly going on, but I suspect this is a poliastro 0.13.1 bug. When it exports data
    # to CZML, the timezone is messed up. I've fixed very similar problem in Poliastro. The patch for this been accepted
    # and will be released in upcoming 0.14.0.
    # This is how it looks like:  "availability": "2020-04-12T09:01:03Z/2020-04-12T09:17:06Z",
    # This is how it SHOULD look: "availability": "2020-04-12T09:01:03/2020-04-12T09:17:06",

    txt = txt.replace("Z/", "/") # replace Z/ with /
    txt = txt.replace('Z"', '"') # replace Z" with "

    txt += "var dataSourcePromise = viewer.dataSources.add(Cesium.CzmlDataSource.load(czml));"

    return txt

def export2cesium(outfile, imgfile, aos, los, aos_list, los_list, methods, tle1, tle2):
    """
    Exports all data to JavaScript that's usable with Cesium.

    Parameters
    ==========
    outfile - output filename to write content to
    imgfile - name of the input PNG file (currently not used yet)
    aos - aquisition of signal, in timedate format
    los - loss of signal, in timedate format
    lla_aos - an array of three coords that specify subsat point in LLA notation (longitude, lattitude, altitude)
    lla_aos - an array of three coords that specify subsat point in LLA notation (longitude, lattitude, altitude)
    tle1 - first line of TLE data
    tle2 - second line of TLE data
    """

    if len(methods) != len(aos_list) or len(methods) != len(los_list):
        raise Exception("Incorrect parameter size: len(methods) != len(aos_list) != len(los_list)")

    txt =  cesium_preamble()
    txt += export2cesium_tle(tle1, tle2, "satname", aos, los)
    txt += "\n\n"

    for i in range(0, len(methods)):
        txt += " // points %d out of %d" % (i, len(methods))
        txt += export2cesium_point(aos_list[i], "AOS:" + str(aos) + ", method " + methods[i])
        txt += export2cesium_point(los_list[i], "LOS:" + str(los) + ", method " + methods[i])

    f = open(outfile, "w")
    f.write(txt)
    f.close()

    print("Georeference data exported to %s" % outfile)

def get_ssp(lla):
    return [ lla[0], lla[1], 0 ]

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

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

    # Ok, we have sat position at time of AOS and LOS. The tricky part here is those are in
    # Earth-Centered Intertial (ECI) reference system. The model used is TEME (True equator,
    # mean equinox). Let's calculate the SSP using three methods:
    # - simple (that assumes spherical Earth)
    # - oblate Earth (uses passed ellipsoid, WGS84 in this case)
    # - using pymap3d lib (which is most precise)

    print("AOS: ECI[x=%f, y=%f, z=%f]" % (pos1[0], pos1[1], pos1[2]))

    aos1 = teme2geodetic_spherical(pos1[0], pos1[1], pos1[2], d1)
    aos2 = teme2geodetic_oblate(pos1[0], pos1[1], pos1[2], d1, ellipsoid_wgs84)
    aos3 = teme2geodetic_pymap3d(pos1[0], pos1[1], pos1[2], d1)

    print("METHOD 1 (spherical Earth): converted to LLA is long=%f lat=%f alt=%f" % (aos1[0], aos1[1], aos1[2]) )
    print("METHOD 2 (oblate Earth):    converted to LLA is long=%f lat=%f alt=%f" % (aos2[0], aos2[1], aos2[2]) )
    print("METHOD 3 (pymap3d):         converted to LLA is long=%f lat=%f alt=%f" % (aos3[0], aos3[1], aos3[2]) )


    los1 = teme2geodetic_spherical(pos2[0], pos2[1], pos2[2], d2)
    los2 = teme2geodetic_oblate(pos2[0], pos2[1], pos2[2], d2, ellipsoid_wgs84)
    los3 = teme2geodetic_pymap3d(pos2[0], pos2[1], pos2[2], d2)

    print("METHOD 1 (spherical Earth): converted to LLA is long=%f lat=%f alt=%f" % (los1[0], los1[1], los1[2]))
    print("METHOD 1 (oblate Earth):    converted to LLA is long=%f lat=%f alt=%f" % (los2[0], los2[1], los2[2]))
    print("METHOD 1 (pymap3d):         converted to LLA is long=%f lat=%f alt=%f" % (los3[0], los3[1], los3[2]))

    # Ok, we have the sat position in LLA format. Getting sub-satellite point is trivial. Just assume altitude is 0.
    aos1 = get_ssp(aos1)
    aos2 = get_ssp(aos2)
    aos3 = get_ssp(aos3)

    los1 = get_ssp(los1)
    los2 = get_ssp(los2)
    los3 = get_ssp(aos3)

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

    # STEP 4: Do we need to detect if it's northbound or southbound fly-over? We can do it easily, by
    # comparing aos_lat (lla1[1]) with los_lat (lla2[1])

    # STEP 5: Export georeferencing data.
    outfile = ".".join(imgname.split('.')[:-1]) + ".js"
    aos_list = [ aos1, aos2, aos3 ]
    los_list = [ los1, los2, los3 ]
    methods  = [ "spherical", "oblate", "pymap3d" ]
    export2cesium(outfile, imgname, d1, d2, aos_list, los_list, methods, tle1, tle2)

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
