from __future__ import division

import math


def normalise_plus_minus_range(value, norm_range):
    """Normalise a value to within a positive and negative range

    **Parameters**

        value
           input value to be limited
        norm_range : {'lat','lon',*number*}
           A number will limit the range between +/- that value. 'lat'
           and 'lon' will limit within +/-90 and 180 respectively.

    """
    norm_range_dict = {'lat': 90, 'lon': 180}
    try:
        valid_range = norm_range_dict[norm_range]
        value %= [-1, 1][value > 0] * valid_range * 2
    except KeyError:
        valid_range = norm_range
    if abs(value) > valid_range:
        return value % ([1, -1][value > 0] * valid_range)
    else:
        return value


class LatLon(object):
    """Implements a latitude Longitude class with conversion to spherical co-ords"""

    def __init__(self, latitude, longitude, unit='degrees'):
        """Performs conversion and stored both degrees and radians.
        The range is limited to within +-90 and +-180 degrees
        """

        if unit in ['deg', 'degrees']:
            _latitude_deg = normalise_plus_minus_range(latitude, 'lat')
            _longitude_deg = normalise_plus_minus_range(longitude, 'lon')
            _latitude_rad = math.radians(_latitude_deg)
            _longitude_rad = math.radians(_longitude_deg)
        elif unit in ['rad', 'radians']:
            _latitude_deg = normalise_plus_minus_range(
                math.degrees(latitude), 'lat')
            _longitude_deg = normalise_plus_minus_range(
                math.degrees(longitude), 'lon')
            _latitude_rad = math.radians(_latitude_deg)
            _longitude_rad = math.radians(_longitude_deg)
        else:
            raise Exception('Incorrect unit specified, please use "deg", "degrees", "rad" or "radians"')
        self._lat = _latitude_deg
        self._lon = _longitude_deg
        self._lat_rad = _latitude_rad
        self._lon_rad = _longitude_rad

        self._default_unit = unit

    @property
    def lat(self):
        """Return Latitude"""
        return self._lat

    @property
    def lon(self):
        """Return Longitude"""
        return self._lon

    @property
    def lon_rad(self):
        """Return Longitude"""
        return self._lon_rad

    @property
    def lat_rad(self):
        """Return Longitude"""
        return self._lat_rad

    def __call__(self):
        return self._lat_lon()

    def __str__(self):
        return "({} N,{} W)".format(*self._lat_lon())

    def _lat_lon(self):
        return self.lat, self.lon

    def convert_spherical(self, altitude, unit='radians'):
        """Converts to spherical co-ordinates

        Returns in radians as standard
        """
        cur_lat_in_radians = self._lat_rad
        equator_radius = 6378137
        flattening = 1 / 298.257223563
        eccentricity_squared = (flattening * (2 - flattening))
        prime_vertical = equator_radius / math.sqrt(1 - eccentricity_squared * math.sin(cur_lat_in_radians)**2)
        p = (prime_vertical + altitude) * math.cos(cur_lat_in_radians)
        z = (prime_vertical * (1 - eccentricity_squared) + altitude) * math.sin(cur_lat_in_radians)
        r = math.sqrt(p**2 + z**2)
        spherical_latitude = math.asin(z / r)
        spherical_longitude = self.lon
        if unit == 'degrees':
            spherical_latitude = math.degrees(spherical_latitude)
            spherical_longitude = math.degrees(spherical_longitude)
        return spherical_latitude, spherical_longitude, r
