from noaa_postproc.imageproc import process_img
import unittest
import pytest

class Georefests(unittest.TestCase):
    def test_ellipsoids(self):

        params = {
            "histogram": True,
            "histogram-adaptive": True,
            "border": True,
            "show": False,
            "write": True,
            "write-left": True,
            "write-right": False,
            "denoise": False,
            "georef": True # Georeference
        }

        process_img("1276.png", params)

        #self.assertIs()
