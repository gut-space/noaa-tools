# This script processes images received from NOAA satellites
#
# Current capabilities:
# - basic sanity checking (if the image is indeed NOAA observation)
# - mark left and right images
# - extract left and right images


import cv2
import sys
from matplotlib import pyplot as plt

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

    # Let's mark borders around the detected images.
    if params['border']:
        img = mark_left(img)
        img = mark_right(img)

    l = extract_left(img)
    r = extract_right(img)

    #show_img(img)
    if params['histogram']:
        l = cv2.equalizeHist(cv2.cvtColor(l, cv2.COLOR_BGR2GRAY))
        r = cv2.equalizeHist(cv2.cvtColor(r, cv2.COLOR_BGR2GRAY))

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

if __name__ == "__main__":
    params = {
        "histogram": True,
        "border": True,
        "show": True,
        "write": False,
        "write-left": True,
        "write-right": True
    }

    process(sys.argv[1], params)
