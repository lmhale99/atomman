# coding: utf-8

# Standard Python libraries
from typing import Optional, Union

# http://www.numpy.org/
import numpy.typing as npt

# atomman imports
from . import IsotropicVolterraDislocation, Stroh, VolterraDislocation
from .. import Box, ElasticConstants

def solve_volterra_dislocation(C: ElasticConstants,
                               burgers: npt.ArrayLike,
                               ξ_uvw: Optional[npt.ArrayLike] = None,
                               slip_hkl: Optional[npt.ArrayLike] = None,
                               transform: Optional[npt.ArrayLike] = None,
                               axes: Optional[npt.ArrayLike] = None,
                               box: Optional[Box] = None,
                               m: Union[str, npt.ArrayLike] = 'x',
                               n: Union[str, npt.ArrayLike] = 'y',
                               cart_axes: bool = False,
                               tol: float = 1e-8) -> VolterraDislocation:
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
    m : str or array-like object, optional
        The 3D Cartesian unit vector to align with the dislocation solution's m-axis,
        i.e. the in-plane direction perpendicular to the dislocation line.  Also
        accepts str values of 'x', 'y', or 'z', in which case the dislocation axis will
        be aligned with the corresponding Cartesian axis.  Default value is 'x'.
    n : str or array-like object, optional
        The 3D Cartesian unit vector to align with the dislocation solution's n-axis,
        i.e. the slip plane normal. Also accepts str values of 'x', 'y', or 'z', in
        which case the dislocation axis will be aligned with the corresponding Cartesian
        axis. Default value is 'y'.
    cart_axes : bool, optional
        Setting this to True will also perform an assertion check that the m- and n-axes
        are both aligned with Cartesian axes. This is a requirement for some of the
        atomic configuration generators. Default value is False as the elastic solution
        by itself does not require the limitation.
    tol : float
        Tolerance parameter used to round off near-zero values.  Default
        value is 1e-8.

    Returns
    -------
    atomman.defect.VolterraDislocation
        The dislocation solution of the appropriate type.
    """

    try:
        return Stroh(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl, transform=transform,
                   axes=axes, box=box, m=m, n=n, cart_axes=cart_axes, tol=tol)
    except ValueError:
        return IsotropicVolterraDislocation(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl,
                                            transform=transform, axes=axes, box=box,
                                            m=m, n=n, cart_axes=cart_axes, tol=tol)
                                            