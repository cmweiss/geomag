'''
Profiling of main class.

Visualisation using SnakeViz
https://jiffyclub.github.io/snakeviz/

'''
import cProfile
import pstats
from datetime import date

from geomag import WorldMagneticModel


def _gen_2d_array(x, y, default=None):
    return [[default] * x for _ in range(y)]


def profile(function):
    cProfile.run(function.__name__ + '()', 'stats')
    p = pstats.Stats('stats')
    p.strip_dirs().sort_stats('time').print_stats(25)


def profile_calc():
    wmm = WorldMagneticModel()
    lat = range(-70, 70)
    lon = range(120)
    world = _gen_2d_array(len(lon), len(lat), date(2015, 1, 1))
    for i, _lon in enumerate(lon):
        for j, _lat in enumerate(lat):
            world[j][i] = wmm.calc_mag_field(_lat, _lon).By
    return world

profile(profile_calc)
