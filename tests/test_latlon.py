import pytest

from geomag.latlon import normalise_plus_minus_range
from geomag.latlon import LatLon


@pytest.mark.parametrize(("lat_lon", "expected_lat_lon"), [
    ((0, 0), (0, 0)),
    ((-91, 181), (89, -179)),
    ((0, 540), (0, 180)),
])
def test_initilisation(lat_lon, expected_lat_lon):
    location = LatLon(*lat_lon)
    assert(location.lat == expected_lat_lon[0])
    assert(location.lon == expected_lat_lon[1])


def test_error_on_wrong_unit():
    with pytest.raises(Exception):
        LatLon((0, 0), (0, 0), 'wrong_unit')


@pytest.mark.parametrize(("value", "norm_range", "expected"), [
    (0, 'lat', 0),
    (0, 'lon', 0),
    (0, 10, 0),
    (20, 10, 0),
    (20, 10, 0),
    (15, 10, -5),
    (91, 'lat', -89),
    (180, 'lon', 180),
    (-180, 'lon', -180),
    (181, 'lon', -179),
    (-181, 'lon', 179),
    (360, 'lon', 0),
    (-360, 'lon', 0),
    (-185, 'lon', 175),
    (540, 'lon', 180),
    (0, 0, 0),
])
def test_norm_range(value, norm_range, expected):
    assert expected == normalise_plus_minus_range(value, norm_range)
