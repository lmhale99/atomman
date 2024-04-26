from fractions import Fraction

import numpy as np

def approx_rational(input, tol=1e-06):
    """
    This finds a rational approximation for all elements in an array.

    Parameters
    ----------
    input : Array-like object
        The array of float values to find rational approximations for.
    tol : float
        The floating point tolerance to use.  The max denominator for the
        approximations will be 1/tol.  Default value is 1e-06.

    Returns
    -------
    numerator : numpy.ndarray
        The identified numerators associated with the approximate rationals
        for each element in input. This will have the same shape as input.
    denominator : numpy.ndarray
        The identified denominators associated with the approximate rationals
        for each element in input. This will have the same shape as input.
    """
    max_denominator = int(1 / tol)

    def find_fraction(input):
        """Finds the approximate rational fraction for a single element"""
        frac = Fraction(input).limit_denominator(max_denominator)
        return frac.numerator, frac.denominator
    vect_find_fraction = np.vectorize(find_fraction, signature='()->(),()')

    numerators, denominators = vect_find_fraction(input)

    return numerators, denominators