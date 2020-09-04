from noaatools import georef
from noaatools.georef import DEG2RAD, RAD2DEG, RE

from math import pi
import unittest
import pytest
from sgp4.api import jday

class Georefests(unittest.TestCase):
    def test_ellipsoids(self):

        x = georef.ellipsoid_wgs84

        self.assertEqual(x.a, 6378.137)
        self.assertEqual(x.finv, 298.257223563)


    def test_julian_date(self):
        # Test set borrowed from here: https://www.astro.umd.edu/~jph/GST_eqn.pdf

        # input: year, month, day, h,m,s, expected output: h,m,s, h(float)
        # Test set borrowed from here: https://www.astro.umd.edu/~jph/GST_eqn.pdf
        expected = [ [2001, 10, 3, 6, 30, 0, 7, 18, 8.32867640, 7.3023135]]

        for case in expected:

            # Convert calendar date to julian date
            jd, fr = jday(case[0], case[1], case[2], case[3], case[4], case[5])

            # Then calculate Greenwich Mean Sidereal Time

            gmst2 = georef.julianDateToGMST2(jd, fr)
            gmst2 = gmst2[0] # julianDateToGSMT2 returns theta and theta_dot, we want just the former
            hms2 = georef.rad2hms(gmst2)

            #print("#### jd=%f fr=%f sum=%f gmst=%f, gmst2=%f hms=%s hms2=%s" % (jd, fr, jd+fr, gmst, gmst2, hms, hms2))

            self.assertAlmostEqual(hms2[0], case[6])
            self.assertAlmostEqual(hms2[1], case[7])
            self.assertAlmostEqual(hms2[2], case[8])
            self.assertAlmostEqual(hms2[3], case[9])

            # TODO: Add checks for julianDateToGMST
            #gmst = georef.julianDateToGMST(jd, fr)
            #hms = georef.rad2hms(gmst)
            #self.assertAlmostEqual(hms[0], case[6])
            #self.assertAlmostEqual(hms[1], case[7])
            #self.assertAlmostEqual(hms[2], case[8])
            #self.assertAlmostEqual(hms[3], case[9])

    def test_rad2hms(self):
        """This test verifies if the RA (expressed in radians can be properly converted to h:m:s)."""

        # Test set borrowed from here: https://www.astro.umd.edu/~jph/GST_eqn.pdf
        expected = [ [7.3023135/24*2*pi, 7, 18, 8.3286, 7.3023135]]

        for case in expected:
            h,m,s,gmst = georef.rad2hms(case[0])

            self.assertAlmostEqual(h, case[1])
            self.assertAlmostEqual(m, case[2])
            self.assertAlmostEqual(s, case[3])
            self.assertAlmostEqual(gmst, case[4])

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
            act2 = georef.julianDateToGMST(case[0], 0.0)

            if type(act) is tuple:
                act = act[0]
            if (abs(exp - act) > 0.000004):
                print("Checking case for jd=%f expected=%f noaatools=%f (diff=%f)" % (case[0], exp, act, abs(exp - act)) )
                self.assertAlmostEqual(exp, act)
            print("JD=%f GMST=%f, GMST2=%f, delta=%f" % (case[0], act2, act, act-act2))

    def test_calc_distance_azimuth(self):
        # lat1, lon1, lat2, lon2 (all in degrees), expected distance in km
        expected = [ [54, 18, 54, 19, 65.43141141], # East
                     [54, 18, 53, 18, 111.31949079], # South
                     [54, 19, 54, 18, 65.43141141], # West
                     [53, 18, 54, 18, 111.31949079], # North
                     [54, 19, 55, 18, 128.72457677],
                     [54, 19, 51, 0, 1325.7042000] ] # Gdansk to London

        print("")

        for c in expected:
            self.assertAlmostEqual(georef.calc_distance(c[0]*DEG2RAD, c[1]*DEG2RAD, c[2]*DEG2RAD, c[3]*DEG2RAD)*RE, c[4])
            self.assertAlmostEqual(georef.distance_apt(c[0]*DEG2RAD, c[1]*DEG2RAD, c[2]*DEG2RAD, c[3]*DEG2RAD)*RE, c[4])

            az1 = float('nan')
            az2 = float('nan')
            az1 = georef.azimuth_apt(c[0]*DEG2RAD, c[1]*DEG2RAD, c[2]*DEG2RAD, c[3]*DEG2RAD)*RAD2DEG
            az2 = georef.calc_azimuth([c[0]*DEG2RAD, c[1]*DEG2RAD], [c[2]*DEG2RAD, c[3]*DEG2RAD])
            try:
                az3 = georef.azimuth_apt(c[0], c[1], c[2], c[3])
            except ValueError:
                print("azimuth_apt failed for the following input: %s" % c)
            try:
                az4 = georef.calc_azimuth([c[0], c[1]], [c[2], c[3]])
            except BadValue:
                print("!!!! calc_azimuth failed for the following input: %s" % c)


            print("%s az1=%f az2=%f\n" % (c, az1, az2))

    def test_longitude_trunc(self):
        self.assertAlmostEqual(georef.longitude_trunc(0), 0)
        self.assertAlmostEqual(georef.longitude_trunc(pi), pi)
        self.assertAlmostEqual(georef.longitude_trunc(-pi), -pi)
        self.assertAlmostEqual(georef.longitude_trunc(2*pi), 0)
        self.assertAlmostEqual(georef.longitude_trunc(2*pi + 0.1), 0.1)
        self.assertAlmostEqual(georef.longitude_trunc(-pi), -pi)
        self.assertAlmostEqual(georef.longitude_trunc(-pi - 1), pi - 1)
        self.assertAlmostEqual(georef.longitude_trunc(100*pi), 0)

