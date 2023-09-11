from enum import Enum
from collections import namedtuple

# Conversion between radians and degrees
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321

RE = 6378.137  # Earth radius (in km)


class Method(Enum):
    SPHERICAL = 1
    OBLATE = 2
    PYMAP3D = 3


# This defines ellipsoid (a = equatorial radius in km, finv = inverse flattening)
Ellipsoid = namedtuple('Ellipsoid', "a finv")

# Source: https://en.wikipedia.org/wiki/Earth_ellipsoid#Historical_Earth_ellipsoids
ellipsoid_wgs84 = Ellipsoid(a=RE, finv=298.257223563)

# Sat image aquisition and image transmission is not instanteneous. There is some delay. It's small, but
# bacause the sat moves at speeds of more than 7km/s, even fractions of second make a big deal of a difference.
# Not sure if this data is available anywhere. I think it'll have to be picked with trial-and-error.
# This is expressed in seconds.
NOAA_PROCESSING_DELAY = 0.0

# This is the FOV (field of view) angle for the observation instrument on NOAA sats for each side.
# In simple words this is how far the camera can look sideways. The total viewing angle is this multipled by 2.
# STK uses 54.3 degrees, however this source claims the FOV to be 55.37 and the swath being 2900km.
# https://directory.eoportal.org/web/eoportal/satellite-missions/n/noaa-poes-series-5th-generation
AVHRR_FOV = 55.37
