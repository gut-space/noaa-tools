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

def show(file: str):
    img = cv2.imread(file)

    ok, error = checkimg(img)
    if not ok:
        print("ERROR: %s" % error)
        return False

    img = mark_left(img)
    img = mark_right(img)
    l = extract_left(img)
    r = extract_right(img)
    #show_img(img)
    #show_img(l)
    #show_img(r)

    fname = file.split('.')
    fname = '.'.join(fname[:-1])

    cv2.imwrite(fname + '-annot.png', img)

    cv2.imwrite(fname + '-left.png', l)
    cv2.imwrite(fname + '-right.png', r)
    print("fname=[%s]" % fname)

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

def show_img(image):
    plt.imshow(image, cmap = 'gray', interpolation = 'bicubic')
    #Uncomment this to hide X, Y axis values
    #plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis

    #plt.plot([200,300,400],[200,200,300],'c', linewidth=5)
    plt.show()

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
    show(sys.argv[1])