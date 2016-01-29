'''
Plots declination for +-70deg latitude and 260deg longditude
'''
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
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

plot_levels = np.arange(-90, 90, 5)

map = Basemap(llcrnrlon=lon[0], llcrnrlat=lat[0], urcrnrlon=lon[-1], urcrnrlat=lat[-1], projection='mill')
map_x, map_y = map(np.asarray(x), np.asarray(y))
contour_set = map.contour(map_x, map_y, world, plot_levels, linewidths=3, cmap=cm.jet)
map.drawcoastlines(linewidth=1.25)
map.fillcontinents(color='0.8')
map.drawparallels(np.arange(lat[0],lat[-1],10),labels=[1,1,0,0])
map.drawmeridians(np.arange(lon[0],lon[-1],20),labels=[0,0,0,1])

plt.clabel(contour_set, fontsize='xx-small', inline_spacing=1.5, fmt='%1.0f')
plt.suptitle('Magnetic Declination (degrees)')
plt.show()
