""" This module contains functions to process NOAA images, mainly drawing grid on them."""

import cv2
# from matplotlib import pyplot as plt

LEFT_BEGIN_COLUMN = 84
LEFT_END_COLUMN = 993

RIGHT_BEGIN_COLUMN = 1124
RIGHT_END_COLUMN = 2033


def process_img(file: str, params):
    """
    Processes PNG image received from NOAA sat.

    Parameters
    ==========
    file: str - filename of the PNG file to be processed
    params : dict - specifies which processing should be done

    Params should have the following format:
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

    Returns
    =======
    True
    """

    img = cv2.imread(file)

    ok, error = checkimg(img)
    if not ok:
        print(f"ERROR: {error}")
        return False

    # img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Denoising (gives poor results so far, need to tweak the parameters)
    if params['denoise']:
        img_d = cv2.fastNlMeansDenoisingColored(img, None, 3, 10, 7, 21)
        show2img(img, img_d)
        img = img_d

    # Let's mark borders around the detected images.
    if params['border']:
        img = mark_left(img)
        img = mark_right(img)

    left = extract_left(img)
    right = extract_right(img)

    # show_img(img)
    if params['histogram']:
        left = cv2.cvtColor(left, cv2.COLOR_BGR2GRAY)
        right = cv2.cvtColor(right, cv2.COLOR_BGR2GRAY)

        if params['histogram-adaptive']:
            clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
            left = clahe.apply(left)
            right = clahe.apply(right)
        else:
            left = cv2.equalizeHist(left)
            right = cv2.equalizeHist(right)

    # Processing done. Let's display them
    if params['show']:
        if params['write-left']:
            show_img(left, title="left", wait=False)
        if params['write-right']:
            show_img(right, title="right", wait=False)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    process_img_write(file, params, img, left, right)

    return True


def process_img_write(file, params, img, left, right):
    """ Writes output files."""

    # Time to write output files.
    fname = file.split('.')
    fname = '.'.join(fname[:-1])

    if params['write']:
        cv2.imwrite(fname + '-annot.png', img)

    if params['write-left']:
        cv2.imwrite(fname + '-left.png', left)

    if params['write-right']:
        cv2.imwrite(fname + '-right.png', right)


def mark_left(img):
    """ Marks left image with a red rectangle."""
    # width and number of channels are ignored.
    height, _, _ = img.shape

    # These are some magic number. The left image starts at pixel 28 and ends on pixel 992.
    return cv2.rectangle(img, (LEFT_BEGIN_COLUMN, 0), (LEFT_END_COLUMN - 1, height - 1), (255, 0, 0))


def mark_right(img):
    """ Marks right image with a green rectangle."""

    height, _, _ = img.shape
    # These are some magic number. The right image starts at pixel 112 and ends on pixel 2033.
    return cv2.rectangle(img, (RIGHT_BEGIN_COLUMN, 0), (RIGHT_END_COLUMN - 1, height - 1), (0, 255, 0))


def show2img(img1, img2, title1="before", title2="after"):
    """
    Shows two images, side by side.

    Parameters
    ==========
    img1, img2 - image one and two
    title1, title2 - title to be displayed
    """
    cv2.imshow(title1, img1)
    cv2.imshow(title2, img2)
    cv2.waitKey(0)
    cv2.destroyAllWindows()


def show_img(img, title="image", wait=True):
    """
    Displays an image. Depending on the wait value, it returns immediately or waits for a key press.
    wait = False is useful when you want to display the first image of several.
    """

    # Let's use simple showing using cv2
    cv2.imshow(title, img)
    if wait:
        cv2.waitKey(0)

    # The alternative is to use pyplot from mathplotlib
    # plt.imshow(image, cmap = 'gray', interpolation = 'bicubic')
    # Uncomment this to hide X, Y axis values
    # plt.xticks([]), plt.yticks([])  # to hide tick values on X and Y axis

    # plt.plot([200,300,400],[200,200,300],'c', linewidth=5)
    # plt.show()


def extract_left(img):
    """ Extracts left image from the full image."""
    _, height, _ = img.shape
    return img[0:height - 1, LEFT_BEGIN_COLUMN:LEFT_END_COLUMN]


def extract_right(img):
    """ Extracts right image from the full image."""
    _, height, _ = img.shape
    return img[0:height - 1, RIGHT_BEGIN_COLUMN:RIGHT_END_COLUMN]


def checkimg(img):
    """ do basic sanity check on the image and see if it looks like NOAA image."""
    height, width, _ = img.shape
    # Check if this looks like NOAA image at all.

    # Check 1: verify it has width of 2080 pixels
    if not height:
        return False, "Image has 0 height."
    if width != 2080:
        return False, "Image has incorrect width: %d, expected 2080 pixels" % width

    return True, ""
