# coding: utf-8

# atomman imports
from . import IsotropicVolterraDislocation, Stroh

def solve_volterra_dislocation(C, burgers, ξ_uvw=None, slip_hkl=None,
                               transform=None, axes=None, box=None,
                               m=[1,0,0], n=[0,1,0], tol=1e-8):
    """
    Wrapper function for generating VolterraDislocation children classes
    that provide linear elastic solutions for straight dislocations.
    
    Parameters
    ----------
    C : atomman.ElasticConstants
        The medium's elastic constants.
    burgers : array-like object
        The dislocation's Burgers vector.
    ξ_uvw : array-like object
        The Miller crystal vector associated with the dislocation's line
        direction.  Must be given with slip_hkl to identify the
        transformation matrix to use on C and burgers.
    slip_hkl : array-like object
        The Miller plane indices associated with the dislocation's slip
        plane.  Must be given with ξ_uvw to identify the
        transformation matrix to use on C and burgers.
    transform : array-like object, optional
        A 3x3 set of orthogonal Cartesian vectors that define the
        transformation matrix to use on C and burgers to convert from the
        standard (unit cell) and dislocation orientations.  The 3 vectors
        will automatically be converted into unit vectors.  Using this is
        an alternative to using ξ_uvw and slip_hkl.
    axes : array-like object, optional
        Same as transform.  Retained for backwards compatibility.
    box : atomman.Box, optional
        The unit cell's box that crystal vectors are taken with respect to.
        If not given, will use a cubic box with a=1 meaning that burgers,
        ξ_uvw and slip_hkl will be interpreted as Cartesian vectors.
    m : array-like object, optional
        The m unit vector for the solution.  m, n, and u (dislocation
        line) should be right-hand orthogonal.  Default value is [1,0,0]
        (x-axis).
    n : array-like object, optional
        The n unit vector for the solution.  m, n, and u (dislocation
        line) should be right-hand orthogonal.  Default value is [0,1,0]
        (y-axis). n is normal to the dislocation slip plane.
    tol : float
        Tolerance parameter used to round off near-zero values.  Default
        value is 1e-8.

    Returns
    atomman.defect.VolterraDislocation
        The dislocation solution of the appropriate type.
    """

    try:
        return Stroh(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl, transform=transform,
                   axes=axes, box=box, m=m, n=n, tol=tol)
    except:
        return IsotropicVolterraDislocation(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl,
                                            transform=transform, axes=axes, box=box,
                                            m=m, n=n, tol=tol)
                                            