'''
Plots declination for +-70deg latitude and 260deg longditude

TODO: Add Basemap example.
'''
import matplotlib.pyplot as plt
from datetime import date
import numpy as np

from geomag import WorldMagneticModel


def _gen_2d_array(size_x, size_y, default=None):
    return [[default] * size_x for _ in range(size_y)]

wmm = WorldMagneticModel()
lat = range(-70, 70)
lon = range(360)
len_lat = len(lat)
len_lon = len(lon)
world = _gen_2d_array(len_lon, len_lat, date(2015, 1, 1))
x = _gen_2d_array(len_lon, len_lat, date(2015, 1, 1))
y = _gen_2d_array(len_lon, len_lat, date(2015, 1, 1))
for i, _lon in enumerate(lon):
    for j, _lat in enumerate(lat):
        world[j][i] = wmm.calc_mag_field(_lat, _lon).declination
        x[j][i] = _lon
        y[j][i] = _lat


plot_levels = np.arange(-90, 90, 10)
fig = plt.figure()
CS = plt.contour(x, y, world, plot_levels)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel('East (degrees)')
plt.ylabel('North (degrees)')
fig.suptitle('Magneteic Declination')
plt.show()
