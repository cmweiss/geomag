from __future__ import division

import math
try:
    import functools32 as functools
except ImportError:
    import functools
import operator


def _gen_2d_array(size_x, size_y, default=None):
    return [[default] * size_x for _ in range(size_y)]


def kronecker_delta(j, i):
    """ Kronecker delta is defined as Iji = 1 if i = j and Iji = 0 otherwise"""
    return 1 if j == i else 0


def double_factorial(n):
    """Returns the double factorial (n!!)

    Input must be an integer or greater than 0 otherwise ValueError is raised

    """
    if int(n) != n:
        raise ValueError('n must be an integer')
    elif n < 0:
        raise ValueError('n must be greater than 0')
    start_number = 2 - (n % 2)
    return functools.reduce(operator.mul, range(start_number, n + 1, 2))


def schmidt_quasi_normalisation(max_n):
    """Returns an array of the Schmidt Quasi-normalised values

    Array is symmetrical about the diagonal
    """
    schmidt = _gen_2d_array(max_n, max_n, 0.0)
    for n in range(max_n):
        for m in range(n + 1):
            if n == 0:
                # This is a bit of a hack to get round 2n-1 evaluating to
                # -1 and erroring when it should be returning 1
                double_fact = 1.0
            else:
                double_fact = double_factorial(2 * n - 1)
            k_delta = kronecker_delta(0, m)
            n_minus_m_fact = math.factorial(n - m)
            n_plus_m_fact = math.factorial(n + m)
            schmidt[m][n] = math.sqrt(
                ((2 - k_delta) * n_minus_m_fact) / n_plus_m_fact) * double_fact / n_minus_m_fact
    return schmidt


def tuplize_array(array):
    return tuple(tuple(line) for line in array)


def recursion_constants(max_n):
    """ Calculates the values useful for performing the scalar potential algorithm """
    k = _gen_2d_array(max_n, max_n, 0)
    for n in range(1, max_n):
        for m in range(n + 1):
            k[m][n] = (((n - 1) * (n - 1)) - (m**2)) / ((2.0 * n - 1) * (2.0 * n - 3.0))
    return tuplize_array(k)


@functools.lru_cache(maxsize=None)
def associated_polynomials(cos_theta, sin_theta, max_n, max_m, k=None):
    """Specific legendre legrende_polynomials for application in geomagnetic calculation

    Recursion Constants (k) can be pre-calculated to improve performance, otherwise the
    function will perform this for you

    """
    if k is None:
        k = recursion_constants(max_n)
    associated_poly = _gen_2d_array(max_n, max_m, 0.0)
    derivative_associated_poly = _gen_2d_array(max_n, max_m, 0.0)
    associated_poly[0][0] = 1
    for n in range(1, max_n):
        for m in range(n + 1):
            if m == n:
                previous_ass_poly = associated_poly[m - 1][n - 1]
                associated_poly[m][n] = cos_theta * previous_ass_poly
                derivative_associated_poly[m][n] = (cos_theta * derivative_associated_poly[m - 1][n - 1]
                                                    + sin_theta * previous_ass_poly)
            else:
                previous_ass_poly = associated_poly[m][n - 1]
                current_k = k[m][n]
                associated_poly[m][n] = sin_theta * previous_ass_poly - current_k * associated_poly[m][n - 2]
                derivative_associated_poly[m][n] = (sin_theta * derivative_associated_poly[m][n - 1]
                                                    - cos_theta * previous_ass_poly
                                                    - current_k * derivative_associated_poly[m][n - 2])
                
    return associated_poly, derivative_associated_poly


def legrende_polynomials(max_n, v):
    """ Returns the legendre polynomials
    Based on the algorithm here:
    http://uk.mathworks.com/help/symbolic/mupad_ref/orthpoly-legendre.html
    """
    poly = [1, v] + [0] * (max_n - 1)
    for n in range(2, max_n + 1):
        poly[n] = (2 * n - 1) / n * v * poly[n - 1] - (n - 1) / n * poly[n - 2]
    return poly


def multiple_angle(phi, max_n):
    """ Returns the values of sin(nx) and cos(nx) series

    """
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    sin_n = [0] * max_n
    cos_n = [1] * max_n
    sin_n[1] = sin_phi
    cos_n[1] = cos_phi
    for i in range(2, max_n):
        last_sin = sin_n[i - 1]
        last_cos = cos_n[i - 1]
        sin_n[i] = sin_phi * last_cos + cos_phi * last_sin
        cos_n[i] = cos_phi * last_cos - sin_phi * last_sin
    return tuple(sin_n), tuple(cos_n)


def scalar_potential(coeff, phi, theta, max_n, max_m, radial_alt, ref_radius=6371200, k=None):
    """

    ref_radius defaults to the radius of the earth in m.
    """
    cos_theta = math.cos(theta)
    if cos_theta is 0:
        cos_theta += 0.00000001
    sin_theta = math.sin(theta)
    sin_n, cos_n = multiple_angle(phi, max_n)
    legendre_poly, legendre_poly_derivative = associated_polynomials(cos_theta, sin_theta, max_n, max_m, k)
    a_over_r_pow = _altitude_ratios(ref_radius, radial_alt, max_n)
    magnetic_field = _calc_scalar_potential(
                                            coeff,
                                            cos_theta,
                                            sin_n,
                                            cos_n,
                                            legendre_poly,
                                            legendre_poly_derivative,
                                            a_over_r_pow,
                                            max_n)
    return magnetic_field


def _calc_scalar_potential(coeff, cos_theta, sin_n, cos_n, legendre_poly, legendre_poly_der, a_over_r_pow, max_n):
    """ Calculates the partial scalar potential values"""
    B_r = 0
    B_theta = 0
    B_phi = 0
    for n in range(1, max_n):
        current_a_over_r_power = a_over_r_pow[n]
        for m in range(n + 1):
            current_cos_n = cos_n[m]
            current_sin_n = sin_n[m]
            g_cos = coeff[m][n] * current_cos_n
            g_sin = coeff[m][n] * current_sin_n
            h_sin = coeff[n][m - 1] * current_sin_n
            h_cos = coeff[n][m - 1] * current_cos_n
            B_r += current_a_over_r_power * (n + 1) * (g_cos + h_sin) * legendre_poly[m][n]
            B_theta -= current_a_over_r_power * (g_cos + h_sin) * legendre_poly_der[m][n]
            B_phi -= current_a_over_r_power * m * (-g_sin + h_cos) * legendre_poly[m][n]
    try:
        B_phi *= 1 / cos_theta
    except ZeroDivisionError:
        B_phi = B_phi
    return B_r, B_theta, B_phi


@functools.lru_cache(maxsize=None)
def _altitude_ratios(ref_radius, altitude_radius, n):
    a_over_r = ref_radius / altitude_radius
    return tuple(a_over_r**(i + 2) for i in range(n))
