# coding: utf-8

# atomman imports
from . import IsotropicVolterraDislocation, Stroh

def solve_volterra_dislocation(C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8, style=None):
    """
    Wrapper function for generating VolterraDislocation children classes
    that provide linear elastic solutions for straight dislocations.
    
    Parameters
    ----------
    C : atomman.ElasticConstants
        The medium's elastic constants.  Must be isotropic.
    burgers : array-like object
        The dislocation's Cartesian Burgers vector.
    axes : array-like object, optional
        3x3 set of rotational axes for the system. If given, C and burgers
        will be transformed using axes.
    m : array-like object, optional
        The m unit vector for the solution.  m, n, and u (dislocation
        line) should be right-hand orthogonal.  Default value is [1,0,0]
        (x-axis).
    n : array-like object, optional
        The n unit vector for the solution.  m, n, and u (dislocation
        line) should be right-hand orthogonal.  Default value is [0,1,0]
        (y-axis). n is normal to the dislocation slip plane.
    tol : float
        Tolerance parameter used to check for compatibility of the other
        parameters.  Default value is 1e-8.
    style : str or None, optional
        Indicates which dislocation solution style to use.  Accepted values are
        'isotropic', 'anisotropic', or None.  If None (default), will select
        isotropic or anisotropic based on the supplied elastic constants.

    Returns
    atomman.defect.VolterraDislocation
        The dislocation solution of the appropriate type.
    """

    if style is None:
        if C.is_normal('isotropic'):
            style = 'isotropic'
        else:
            style = 'anisotropic'
    if style == 'isotropic':
        return IsotropicVolterraDislocation(C, burgers, axes=axes, m=m, n=n, tol=tol)
    elif style == 'anisotropic':
        return Stroh(C, burgers, axes=axes, m=m, n=n, tol=tol)
    else:
        raise ValueError('Invalid style: must be "isotropic" or "anisotropic"')

    
