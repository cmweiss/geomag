"""
by Christopher Weiss cmweiss@gmail.com
and Todd Dembrey todd.dembrey@gmail.com

Adapted from the geomagc software and World Magnetic Model of the NOAA
Satellite and Information Service, National Geophysical Data Center
http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml

Suggestions for improvements are appreciated.

"""
from .world_magnetic_model import WorldMagneticModel

__all__ = ['WorldMagneticModel']


# The following functions are for backwards compatability

def declination(*args, **kargs):
    return WorldMagneticModel().calc_mag_field(*args, **kargs).declination


def heading(hdg, *args, **kargs):
    return WorldMagneticModel().calc_mag_field(*args, **kargs).mag_heading(hdg)
