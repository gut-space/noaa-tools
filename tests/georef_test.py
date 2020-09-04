from noaatools import georef
from math import pi
import unittest
import pytest

class Georefests(unittest.TestCase):
    def test_ellipsoids(self):

        x = georef.ellipsoid_wgs84

        self.assertEqual(x.a, 6378.137)
        self.assertEqual(x.finv, 298.257223563)


    def test_theta(self):
        """This test verifies if julian date conversion to Greenwich Mean Sidreal Time is correct."""

        # Data generated using https://github.com/skyfielders/python-skyfield/blob/master/skyfield/sgp4lib.py
        # for x in np.linspace(2458978.5, 2458980.5, 20):
        #     print("[%f,%f]," % (x, sgp4lib.theta_GMST1982(x)[0]))
        # 2458978.5 is 8 May 2020, 2458980.500000 is 11 May 2020
        expected = [ [2458978.500000,3.966616],
                    [2458978.605263,4.629814],
                    [2458978.710526,5.293013],
                    [2458978.815789,5.956212],
                    [2458978.921053,0.336225],
                    [2458979.026316,0.999424],
                    [2458979.131579,1.662623],
                    [2458979.236842,2.325822],
                    [2458979.342105,2.989020],
                    [2458979.447368,3.652219],
                    [2458979.552632,4.315418],
                    [2458979.657895,4.978616],
                    [2458979.763158,5.641815],
                    [2458979.868421,0.021828],
                    [2458979.973684,0.685027],
                    [2458980.078947,1.348226],
                    [2458980.184211,2.011425],
                    [2458980.289474,2.674624],
                    [2458980.394737,3.337822],
                    [2458980.500000,4.001021] ]

        for case in expected:
            exp = case[1]
            act = georef.julianDateToGMST2(case[0], 0.0)

            if type(act) is tuple:
                act = act[0]
            if (abs(exp - act) > 0.000004):
                print("Checking case for jd=%f expected=%f noaatools=%f (diff=%f)" % (case[0], exp, act, abs(exp - act)) )
                self.assertAlmostEqual(exp, act)

    def test_calc_distance(self):
        self.assertAlmostEqual(georef.calc_distance(54, 18, 54, 19), 65.43141141)
        self.assertAlmostEqual(georef.calc_distance(54, 18, 54, 17), 65.43141141)
        self.assertAlmostEqual(georef.calc_distance(54, 18, 53, 18), 111.31949079)
        self.assertAlmostEqual(georef.calc_distance(54, 19, 55, 18), 128.72457677)
        self.assertAlmostEqual(georef.calc_distance(54, 19, 51, 0), 1325.7042000) # Gdansk to London


    def test_longitude_trunc(self):
        self.assertAlmostEqual(georef.longitude_trunc(0), 0)
        self.assertAlmostEqual(georef.longitude_trunc(pi), pi)
        self.assertAlmostEqual(georef.longitude_trunc(-pi), -pi)
        self.assertAlmostEqual(georef.longitude_trunc(2*pi), 0)
        self.assertAlmostEqual(georef.longitude_trunc(2*pi + 0.1), 0.1)
        self.assertAlmostEqual(georef.longitude_trunc(-pi), -pi)
        self.assertAlmostEqual(georef.longitude_trunc(-pi - 1), pi - 1)
        self.assertAlmostEqual(georef.longitude_trunc(100*pi), 0)

