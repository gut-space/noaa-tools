from noaatools.georef import usage, ellipsoid_wgs84
import unittest
import pytest

class Georefests(unittest.TestCase):
    def test_ellipsoids(self):

        x = ellipsoid_wgs84

        self.assertEqual(x.a, 6378.137)
        self.assertEqual(x.finv, 298.257223563)

