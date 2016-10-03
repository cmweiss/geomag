# -*- coding: utf-8 -*-
"""
© 2012 - 2016 Xample Sàrl

Author: Stephane Poss
Date: 03.10.16
"""
import unittest
import datetime

from geomag import GeoMag


class GeoMagTest(unittest.TestCase):
    def setUp(self):
        d1 = datetime.date(2015, 1, 1)
        d2 = datetime.date(2017, 7, 2)

        self.test_values = (
            # date, alt, lat, lon, var
            (d1, 0, 80, 0, -3.85),
            (d1, 0, 0, 120, 0.57),
            (d1, 0, -80, 240, 69.81),
            (d1, 328083.99, 80, 0, -4.27),
            (d1, 328083.99, 0, 120, 0.56),
            (d1, 328083.99, -80, 240, 69.22),
            (d2, 0, 80, 0, -2.75),
            (d2, 0, 0, 120, 0.32),
            (d2, 0, -80, 240, 69.58),
            (d2, 328083.99, 80, 0, -3.17),
            (d2, 328083.99, 0, 120, 0.32),
            (d2, 328083.99, -80, 240, 69.00),
        )

    def test_declination(self):
        gm = GeoMag()
        for values in self.test_values:
            calcval = gm.GeoMag(values[2], values[3], values[1], values[0])
            self.assertAlmostEqual(values[4], calcval.dec, 2, 'Expected %s, result %s' % (values[4], calcval.dec))
