import argparse
from noaatools import export_csv
from noaatools import export_js
from noaatools import georef
from datetime import timedelta

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


def usage():
    print(USAGE)


def method2enum(str):
    if str == "SPHERICAL":
        return georef.Method.SPHERICAL
    if str == "OBLATE":
        return georef.Method.OBLATE
    if str == "PYMAP3D":
        return georef.Method.PYMAP3D
    raise Exception("ERROR: Unable to understand %s, allowed values are SPHERICAL, OBLATE, PYMAP3D" % str)


if __name__ == "__main__":

    # These are the default parameters. Those are taken from satnogs project, observation 1276.
    tle1 = '1 28654U 05018A   20098.54037539  .00000075  00000-0  65128-4 0  9992'
    tle2 = '2 28654  99.0522 154.2797 0015184  73.2195 287.0641 14.12501077766909'

    aos = '2020-04-12Z09:01:03.063476'
    los = '2020-04-12Z09:17:06.466954'

    satname = 'NOAA 18'
    imgname = 'data/1276.png'

    method = "PYMAP3D"

    parser = argparse.ArgumentParser(APP_NAME)
    parser.add_argument("--aos", help="Specify Acquision of Signal as timestamp (default: %(default)s)", default=aos, type=str)
    parser.add_argument("--los", help="Specify Loss of signal as timestamp (default: %(default)s)", default=los, type=str)
    parser.add_argument("--tle1", help="First line of TLE (default: %(default)s)", default=tle1, type=str)
    parser.add_argument("--tle2", help="Second line of TLE (default: %(default)s)", default=tle2, type=str)
    parser.add_argument(
        "--method",
        help="Specify orbital position calculation method (SPHERICAL, OBLATE, PYMAP3D, default: %(default)s)",
        default=method,
        type=str)
    parser.add_argument("--satname", help="Specifies satellite name (default: %(default)s)", default=satname, type=str)
    parser.add_argument("--file", help="Specifies output file pattern (default: %(default)s)", default=imgname, type=str)

    args = parser.parse_args()
    aos = args.aos
    los = args.los
    tle1 = args.tle1
    tle2 = args.tle2
    method = method2enum(args.method)
    satname = args.satname
    imgname = args.file

    # Calculate orbital positions (AOS, LOS), the camera geometry and get the image corners
    d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr = georef.georef(method, tle1, tle2, aos, los)

    # Now export the data to Cesium JavaScript
    outfile_js = ".".join(imgname.split('.')[:-1]) + ".js"
    outfile_csv = ".".join(imgname.split('.')[:-1]) + ".txt"

    # BUG: For some reason the mean anomaly calculated in CZML exporter is off by couple degrees that's roughly equivalent to 5 minutes time.
    # Let's shift time by 5 minutes.
    delta = timedelta(minutes=5)
    d1_czml = d1 - delta
    d2_czml = d2 - delta

    export_js.export2cesium(outfile_js, satname, d1_czml, d2_czml, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, method.name)
    export_csv.export2csv(outfile_csv, satname, d1, d2, aos_lla, los_lla, corner_ul, corner_ur, corner_ll, corner_lr, tle1, tle2, method.name)

    # STEP 6: (possibly outside of this script):
    # - use GDAL library to georeference image (https://pcjericks.github.io/py-gdalogr-cookbook/)
    # - import into QGIS (and follow this tutorial: https://www.qgistutorials.com/en/docs/georeferencing_basics.html)
    # - display the image in Cesium (https://www.cesium.com/docs/cesiumjs-ref-doc/SingleTileImageryProvider.html)
