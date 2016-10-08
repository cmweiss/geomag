from __future__ import division

import math
import os
from datetime import date

from .scalar_potential import scalar_potential, schmidt_quasi_normalisation, recursion_constants
from .latlon import LatLon


def _gen_square_array(size_x, default=None):
    """Creates a square array with x by x elements"""
    return [[default] * size_x for _ in range(size_x)]


def _calculate_decimal_year(date_of_year):
    """
    .total_seconds() call makes the function 2/3 compliant. timedelta division was added in 3.2.
    This does limit the module to 2.7 as the minimum supported
    """
    year = date_of_year.year
    start_of_this_year = date(year, 1, 1)
    days_this_year = (date(year + 1, 1, 1) - start_of_this_year).total_seconds()
    days_into_year = (date_of_year - start_of_this_year).total_seconds()
    return year + days_into_year / days_this_year


def _convert_to_km(value, unit):
    conversion_factor = {'ft': 3280.8399, 'm': 1000, 'km': 1}
    try:
        value_in_km = value / conversion_factor[unit]
    except KeyError:
        raise KeyError('Unknown unit: {unit}')
    return value_in_km

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), 'model_data/WMM.COF')


class MagneticModelData(object):
    """ Class to hold the data and calculations from the file."""
    _last_calculated_datetime = None
    max_order = degree_of_expansion = 12
    array_size = max_order + 1
    coeff = {}
    coefficient = _gen_square_array(array_size, 0.0)
    coefficient_dot = _gen_square_array(array_size, 0.0)
    time_adjusted_coefficients_cache = _gen_square_array(array_size, 0.0)
    
    def __init__(self, file):
        with open(file) as world_magnetic_model_file:
            for line in world_magnetic_model_file:
                linevals = line.strip().split()
                if len(linevals) == 3:
                    self.epoch = float(linevals[0])
                    self.model = linevals[1]
                    self.modeldate = linevals[2]
                elif len(linevals) == 6:
                    degree_n = int(float(linevals[0]))
                    order_m = int(float(linevals[1]))
                    gauss_g = float(linevals[2])
                    gauss_h = float(linevals[3])
                    gauss_g_dot = float(linevals[4])
                    gauss_h_dot = float(linevals[5])
                    self.coefficient[order_m][degree_n] = gauss_g
                    self.coefficient_dot[order_m][degree_n] = gauss_g_dot
                    if order_m != 0:
                        self.coefficient[degree_n][order_m - 1] = gauss_h
                        self.coefficient_dot[degree_n][order_m - 1] = gauss_h_dot
        self._unnormalise_gauss_coefficients()
        
    def _unnormalise_gauss_coefficients(self):
        """ Convert Schmidt normalized Gauss coefficients to unnormalized """
        schmidt_norm = schmidt_quasi_normalisation(self.array_size)
        for n in range(self.array_size):
            for m in range(self.array_size):
                if m <= n:
                    self.coefficient[m][n] = schmidt_norm[m][n] * self.coefficient[m][n]
                    self.coefficient_dot[m][n] = schmidt_norm[m][n] * self.coefficient_dot[m][n]
                else:
                    self.coefficient[m][n] = schmidt_norm[n + 1][m] * self.coefficient[m][n]
                    self.coefficient_dot[m][n] = schmidt_norm[n + 1][m] * self.coefficient_dot[m][n]

    def time_adjust_gauss(self, time):
        """Time adjust the Gauss Coefficients
        
            There is a very basic cache happening here where if the time delta 
            hasn't changed then the previous calculation is returned
        """
        current_delta_time = _calculate_decimal_year(time) - self.epoch
        if self._last_calculated_datetime != current_delta_time:
            self._update_time_coefficients(current_delta_time)
            self._last_calculated_datetime = current_delta_time
        return self.time_adjusted_coefficients_cache

    def _update_time_coefficients(self, delta_time):
        for n in range(1, self.max_order + 1):
            for m in range(0, n + 1):
                self.time_adjusted_coefficients_cache[m][n] = (self.coefficient[m][n] +
                                                               delta_time * self.coefficient_dot[m][n])
                if m != 0:
                    self.time_adjusted_coefficients_cache[n][m - 1] = (self.coefficient[n][m - 1] +
                                                                       delta_time * self.coefficient_dot[n][m - 1])


class WorldMagneticModel(object):
    """Class for calculating geomagnetic variation according to the world magnetic model

    Example Usage:

    >>> from geomag import WorldMagneticModel
    >>> wmm = WorldMagneticModel()
    >>> wmm.calc_mag_field(80,0).declination

    """
    
    radius_earth = 6371200

    _northerly_intensity = None
    _easterly_intensity = None
    _vertical_intensity = None
    _horizontal_intensity = None
    _total_intensity = None
    _declination = None
    _inclination = None
    _grid_variation = None
        
    def __init__(self, world_magnetic_model_filename=DEFAULT_PATH):
        """__init__(self,world_magnetic_model_filename='WMM.COF')
        Loads a file containing the constants for the magnetic model.

        The coefficients for the model are included for the year 2015 if no
        variable is provided.
        File should be in the format of:

        =========    ========    ========    ========    ========    ========
        degree(n)    order(m)    g |mnt0|    h |mnt0|    g |mnt0|    h |mnt0|
        =========    ========    ========    ========    ========    ========

        .. |mnt0| replace:: \ :sup:`m`:sub:`n`\ (t\ :sub:`0`\ )

        """
        self.data = MagneticModelData(world_magnetic_model_filename)
        self.k = recursion_constants(self.data.array_size)
        
    @property
    def grid_variation(self):
        return self._grid_variation

    @property
    def dip(self):
        return self._inclination

    @property
    def declination(self):
        return self._declination

    @property
    def total_intensity(self):
        return self._total_intensity

    @property
    def Bh(self):
        return self._horizontal_intensity

    @property
    def Bx(self):
        return self._northerly_intensity

    @property
    def By(self):
        return self._easterly_intensity

    @property
    def Bz(self):
        return self._vertical_intensity

    def calc_mag_field(self, dlat, dlon, altitude=0, date=date.today(), unit='ft'):
        """calc_mag_field(self, dlat, dlon, altitude=0, date=date.today(), unit='ft')

        Calculates the magnetic field for a given latitude and longitude in decimal degrees.

        **Parameters**

            dlat
               Latitude in degrees
            dlon
               Longitude in degrees
            altitude : optional
               Altitude at which to evaluate magnetic field
            date : datetime.date - optional
               Time will default to today
            unit : {'ft','m','km'} - optional
               Unit for altitude

        """
        altitude_in_km = _convert_to_km(altitude, unit)

        if not -1 <= altitude_in_km <= 850:
            raise ValueError('World Magnetic Model is not valid outside the rage -1 to 850km')

        lat_lon = LatLon(dlat, dlon)
        spherical_latitude, _, radial = lat_lon.convert_spherical(altitude_in_km * 1000)

        b_radius, b_theta, b_phi = scalar_potential(self.data.time_adjust_gauss(date),
                                                    lat_lon.lon_rad,
                                                    spherical_latitude,
                                                    self.data.array_size, 
                                                    self.data.array_size,
                                                    radial, 
                                                    k=self.k)
        # Matching the method in the document
        b_radius = -b_radius
        b_theta = -b_theta
        # /*
        # ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
        # GEODETIC COORDINATES
        # */
        delta_latitude_radians = spherical_latitude - lat_lon.lat_rad
        sin_delta_latitude = math.sin(delta_latitude_radians)
        cos_delta_latitude = math.cos(delta_latitude_radians)

        self._northerly_intensity = b_theta * cos_delta_latitude - b_radius * sin_delta_latitude
        self._easterly_intensity = b_phi
        self._vertical_intensity = b_theta * sin_delta_latitude + b_radius * cos_delta_latitude
        self._horizontal_intensity = math.sqrt(self._northerly_intensity ** 2 + self._easterly_intensity**2)
        self._total_intensity = math.sqrt(self._horizontal_intensity**2 + self._vertical_intensity**2)
        self._declination = math.degrees(math.atan2(self._easterly_intensity, self._northerly_intensity))
        self._inclination = math.degrees(math.atan2(self._vertical_intensity, self._horizontal_intensity))
        self._grid_variation = self._calculate_grid_variation(lat_lon)
        return self

    def _calculate_grid_variation(self, lat_lon):
        """Calculate the magnetic grid variation

        Compute magnetic grid variation if the current
        geodetic position is in the arctic or antarctic
        (i.e. glat > +55 degrees or glat < -55 degrees)
        Otherwise, set magnetic grid variation to 0
        """
        cur_lat = lat_lon.lat
        grid_variation = self._declination
        if cur_lat > 55:
            grid_variation -= lat_lon.lon
        elif cur_lat < -55:
            grid_variation += lat_lon.lon
        grid_variation %= 360
        return grid_variation

    def mag_heading(self, hdg):
        """Calculates the magnetic heading from a true heading.
        hdg = true heading in degrees
        Example Usage:

        >>> wmm.calc_mag_field(0,80).mag_heading(20)

        """
        return (hdg - self.declination + 360.0) % 360

    def field_vectors(self):
        """ Returns the main magnetic components: X, Y and Z"""
        return self.return_array('X', 'Y', 'Z')

    def return_array(self, *variable_list):
        """Populates a return array with the requested variables, where:

        | X = Northerly intensity
        | Y = Easterly intensity
        | Z = Vertical intensity
        | H = Horizontal Intensity
        | F = Total Intensity
        | I = Inclination
        | D = Dip
        | GV = Grid Variation

        Example Usage:

        >>> wmm.calc_mag_field(0,80).return_array('GV','I')

        """
        variable_dict = {
            'X': self.Bx, 'Y': self.By, 'Z': self.Bz, 'H': self.Bh,
            'F': self.total_intensity, 'I': self.dip, 'D': self.dip,
            'GV': self.grid_variation}
        return_list = []
        for variable in list(variable_list):
            return_list.append(variable_dict[variable])
        return return_list
