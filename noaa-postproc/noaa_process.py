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

import cv2
import sys
from datetime import datetime, timezone
from matplotlib import pyplot as plt

from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.api import jday, Satrec

LEFT_BEGIN_COLUMN = 84
LEFT_END_COLUMN = 993

RIGHT_BEGIN_COLUMN = 1124
RIGHT_END_COLUMN = 2033

def process(file: str, params):
    img = cv2.imread(file)

    ok, error = checkimg(img)
    if not ok:
        print("ERROR: %s" % error)
        return False

    #img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Denoising (gives poor results so far, need to tweak the parameters)
    if params['denoise']:
        img_d = cv2.fastNlMeansDenoisingColored(img, None, 3, 10, 7, 21)
        show2img(img, img_d)
        img = img_d

    # Let's mark borders around the detected images.
    if params['border']:
        img = mark_left(img)
        img = mark_right(img)

    l = extract_left(img)
    r = extract_right(img)

    #show_img(img)
    if params['histogram']:
        l = cv2.cvtColor(l, cv2.COLOR_BGR2GRAY)
        r = cv2.cvtColor(r, cv2.COLOR_BGR2GRAY)

        if params['histogram-adaptive']:
            clahe = cv2.createCLAHE(clipLimit = 2.0, tileGridSize=(8,8))
            l = clahe.apply(l)
            r = clahe.apply(r)
        else:
            l = cv2.equalizeHist(l)
            r = cv2.equalizeHist(r)

    # Processing done. Let's display them
    if params['show']:
        if params['write-left']:
            show_img(l, title="left", wait=False)
        if params['write-right']:
            show_img(r, title="right", wait=False)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    # Time to write output files.
    fname = file.split('.')
    fname = '.'.join(fname[:-1])

    if params['write']:
        cv2.imwrite(fname + '-annot.png', img)

    if params['write-left']:
        cv2.imwrite(fname + '-left.png', l)

    if params['write-right']:
        cv2.imwrite(fname + '-right.png', r)

    return True

def mark_left(img):
    # width and number of channels are ignored.
    height, _, _ = img.shape

    print("#height=%s" % height)
    # These are some magic number. The left image starts at pixel 28 and ends on pixel 992.
    return cv2.rectangle(img, (LEFT_BEGIN_COLUMN, 0), (LEFT_END_COLUMN - 1, height - 1), (255,0,0) )

def mark_right(img):

    height, _, _ = img.shape
    # These are some magic number. The right image starts at pixel 112 and ends on pixel 2033.
    return cv2.rectangle(img, (RIGHT_BEGIN_COLUMN, 0), (RIGHT_END_COLUMN - 1, height - 1), (0,255,0))

def show2img(img1, img2, title1 = "before", title2 = "after"):
    cv2.imshow(title1, img1)
    cv2.imshow(title2, img2)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def show_img(img, title="image", wait = True):

    # Let's use simple showing using cv2
    cv2.imshow(title, img)
    if wait:
        cv2.waitKey(0)

    # The alternative is to use pyplot from mathplotlib
    #plt.imshow(image, cmap = 'gray', interpolation = 'bicubic')
    #Uncomment this to hide X, Y axis values
    #plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis

    #plt.plot([200,300,400],[200,200,300],'c', linewidth=5)
    #plt.show()

def extract_left(img):
    _, height, _ = img.shape
    return img[0:height-1, LEFT_BEGIN_COLUMN:LEFT_END_COLUMN]

def extract_right(img):
    _, height, _ = img.shape
    return img[0:height-1, RIGHT_BEGIN_COLUMN:RIGHT_END_COLUMN]

def checkimg(img):
    height, width, _ = img.shape
    # Check if this looks like NOAA image at all.

    # Check 1: verify it has width of 2080 pixels
    if not height:
        return False, "Image has 0 height."
    if width != 2080:
        return False, "Image has incorrect width: %d, expected 2080 pixels" % width

    return True, ""

def georef_naive(tle1, tle2, aos, los):
    """ This is a naive georeferencing method:
        - calculates the sat location at AOS and LOS points (using )
    then calculates distance between them. """

    d1 = datetime.fromisoformat(aos).replace(tzinfo=timezone.utc)
    d2 = datetime.fromisoformat(los).replace(tzinfo=timezone.utc)

    # STEP 1: Calculate sat location at AOS and LOS

    # This approach uses old API 1.x
    sat_old = twoline2rv(tle1, tle2, wgs72)
    pos1_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    pos2_old, _ = sat_old.propagate(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)

    # This approach uses new API 2.x which gives a slightly different results.
    # In case of NOAA, the position is off by less than milimeter
    sat = Satrec.twoline2rv(tle1, tle2)
    jd1, fr1 = jday(d1.year, d1.month, d1.day, d1.hour, d1.minute, d1.second)
    jd2, fr2 = jday(d2.year, d2.month, d2.day, d2.hour, d2.minute, d2.second)
    _, pos1, _ = sat.sgp4(jd1, fr1) # returns error, position and velocity - we care about position only
    _, pos2, _ = sat.sgp4(jd2, fr2)

    # STEP 2: Calculate sub-satellite point at AOS, LOS times
    # T.S. Kelso saves the day *again*: see here: https://celestrak.com/columns/v02n03/

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
