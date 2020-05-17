import sys
import os

APP_NAME = "noaatools"

# This script processes images received from NOAA satellites
USAGE = '''
#
# noaa-cli script
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

from datetime import timedelta
from noaatools import georef
from noaatools import export_js

import argparse

def usage():
    print(USAGE)

if __name__ == "__main__":

    # Let's ignore input parameters and pretend we were asked to georeference observation #1276.
    if len(sys.argv) == 1:
        usage()
        sys.exit(-1)

    parser = argparse.ArgumentParser(APP_NAME)

    tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
    tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

    aos = '2020-04-12Z09:01:03.063476'
    los = '2020-04-12Z09:17:06.466954'

    satname = 'NOAA 18'
    imgname = 'data/1276.png'

    method = georef.Method.SPHERICAL

    # Calculate orbital positions (AOS, LOS), the camera geometry and get the image corners
    d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr = georef.georef(method, tle1, tle2, aos, los)

    # Now export the data to Cesium JavaScript
    outfile = ".".join(imgname.split('.')[:-1]) + ".js"

    # BUG: For some reason the mean anomaly calculated in CZML exporter is off by couple degrees that's roughly equivalent to 5 minutes time.
    # Let's shift time by 5 minutes.
    delta = timedelta(minutes=5)
    d1 -= delta
    d2 -= delta

    export_js.export2cesium(outfile, satname, d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, method.name)

    # STEP 6: (possibly outside of this script):
    # - use GDAL library to georeference image (https://pcjericks.github.io/py-gdalogr-cookbook/)
    # - import into QGIS (and follow this tutorial: https://www.qgistutorials.com/en/docs/georeferencing_basics.html)
    # - display the image in Cesium (https://www.cesium.com/docs/cesiumjs-ref-doc/SingleTileImageryProvider.html)
