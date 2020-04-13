# This script processes images received from NOAA satellites
#
# Current capabilities:
# - basic sanity checking (if the image is indeed NOAA observation)
# - mark left and right images
# - extract left and right images


import cv2
import sys
from matplotlib import pyplot as plt

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
    show_img(img)
    show_img(l)
    show_img(r)
    return True

def mark_left(img):
    # width and number of channels are ignored.
    _, height, _ = img.shape

    # These are some magic number. The left image starts at pixel 28 and ends on pixel 992.
    return cv2.rectangle(img, (84, 0), (992, height), (255,0,0) )

def mark_right(img):

    _, height, _ = img.shape
    # These are some magic number. The right image starts at pixel 112 and ends on pixel 2033.
    return cv2.rectangle(img, (1124, 0), (2033, height), (0,255,0))

    return img

def show_img(image):
    plt.imshow(image, cmap = 'gray', interpolation = 'bicubic')
    #Uncomment this to hide X, Y axis values
    #plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis

    #plt.plot([200,300,400],[200,200,300],'c', linewidth=5)
    plt.show()

def extract_left(img):
    _, height, _ = img.shape
    return img[0:height, 84:992]

def extract_right(img):
    _, height, _ = img.shape
    return img[0:height, 1124:2033]

def checkimg(img):
    height, width, channels = img.shape
    # Check if this looks like NOAA image at all.

    # Check 1: verify it has width of 2080 pixels
    if not height:
        return False, "Image has 0 height."
    if width != 2080:
        return False, "Image has incorrect width: %d, expected 2080 pixels" % width


    return True, ""

if __name__ == "__main__":
    show(sys.argv[1])