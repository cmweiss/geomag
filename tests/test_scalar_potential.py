import pytest

from geomag.scalar_potential import legrende_polynomials, associated_polynomials, double_factorial


@pytest.mark.parametrize(("n", "expected"), [
    (1, 1),
    (2, 2),
    (3, 3),
    (4, 8),
    (5, 15),
])
def test_double_factorial(n, expected):
    assert double_factorial(n) == expected


@pytest.mark.parametrize(("in_array", "expected_poly"), [
    ((1, 1), [1, 1]),
    ((2, 1), [1, 1, 1]),
    ((3, 1), [1, 1, 1, 1]),
    ((1, 2), [1, 2]),
    ((2, 2), [1, 2, 5.5]),
    ((3, 2), [1, 2, 5.5, 17]),
])
def test_polynomials(in_array, expected_poly):
    assert legrende_polynomials(*in_array) == pytest.approx(expected_poly)


@pytest.mark.parametrize(("in_array", "expected_poly"), [
    (({'cos_theta': 1, 'sin_theta': 1, 'max_n': 1, 'max_m': 1}), ([[1]], [[0]])),
])
def test_associated_polynomials(in_array, expected_poly):
    assert associated_polynomials(**in_array) == expected_poly
