import pytest
from collections import namedtuple
from datetime import date

from geomag import WorldMagneticModel

DATE_IN_2015 = date(2015, 1, 1)
DATE_IN_2017 = date(2017, 7, 2)

WMM = WorldMagneticModel()

RESULTS = namedtuple('RESULTS', ['X', 'Y', 'Z', 'H', 'F', 'I', 'D', 'GV'])
CASES = [([80, 0, 0, DATE_IN_2015],
          RESULTS(6627.1, -445.9, 54432.3, 6642.1, 54836, 83.04, -3.85, -3.85)),
         ([0, 120, 0, DATE_IN_2015],
          RESULTS(39518.2, 392.9, -11252.4, 39520.2, 41090.9, -15.89, 0.57, 0.57)),
         ([-80, 240, 0, DATE_IN_2015],
          RESULTS(5797.3, 15761.1, -52919.1, 16793.5, 55519.8, -72.39, 69.81, 309.81)),
         ([80, 0, 328083.99, DATE_IN_2015],
          RESULTS(6314.3, -471.6, 52269.8, 6331.9, 52652, 83.09, -4.27, -4.27)),
         ([0, 120, 328083.99, DATE_IN_2015],
          RESULTS(37535.6, 364.4, -10773.4, 37537.3, 39052.7, -16.01, 0.56, 0.56)),
         ([-80, 240, 328083.99, DATE_IN_2015],
          RESULTS(5613.1, 14791.5, -50378.6, 15820.7, 52804.4, -72.57, 69.22, 309.22)),
         ([80, 0, 0, DATE_IN_2017],
          RESULTS(6599.4, -317.1, 54459.2, 6607, 54858.5, 83.08, -2.75, -2.75)),
         ([0, 120, 0, DATE_IN_2017],
          RESULTS(39571.4, 222.5, -11030.1, 39572, 41080.5, -15.57, 0.32, 0.32)),
         ([-80, 240, 0, DATE_IN_2017],
          RESULTS(5873.8, 15781.4, -52687.9, 16839.1, 55313.4, -72.28, 69.58, 309.58)),
         ([80, 0, 328083.99, DATE_IN_2017],
          RESULTS(6290.5, -348.5, 52292.7, 6300.1, 52670.9, 83.13, -3.17, -3.17)),
         ([0, 120, 328083.99, DATE_IN_2017],
          RESULTS(37585.5, 209.5, -10564.2, 37586.1, 39042.5, -15.7, 0.32, 0.32)),
         ([-80, 240, 328083.99, DATE_IN_2017],
          RESULTS(5683.5, 14808.8, -50163, 15862, 52611.1, -72.45, 69, 309)),
         ]


@pytest.fixture(scope="class", params=CASES)
def setup_class(request):
    print(request)
    request.cls.results = WMM.calc_mag_field(*request.param[0])
    request.cls.expected = request.param[1]


@pytest.mark.usefixtures("setup_class")
class TestGeoMag:

    def test_northerly_intensity(self):
        assert(self.expected.X - self.results.Bx < 1)

    def test_easterly_intensity(self):
        assert(self.expected.Y - self.results.By < 1)

    def test_down_intensity(self):
        assert(self.expected.Z - self.results.Bz < 1)

    def test_total_intensity(self):
        assert(self.expected.F - self.results.total_intensity < 1)

    def test_inclination(self):
        assert(self.expected.I - self.results.dip < 0.01)

    def test_declination(self):
        assert(self.expected.D - self.results.declination < 0.01)

    def test_grid_variation(self):
        assert(self.expected.GV - self.results.grid_variation < 0.01)
