# coding: utf-8
# Standard Python libraries
from typing import Optional

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from .. import Box
from ..tools import miller

def dislocation_system_transform(ξ_uvw: npt.ArrayLike,
                                 slip_hkl: npt.ArrayLike,
                                 m: Optional[npt.ArrayLike] = None,
                                 n: Optional[npt.ArrayLike] = None,
                                 box: Optional[Box] = None,
                                 tol: float = 1e-8) -> np.ndarray:
    """
    Utility function for getting the transformation matrix to the dislocation system
    based on crystal orientation vectors.
    
    Parameters
    ----------
    ξ_uvw : array-like object
        The Miller crystal vector associated with the dislocation's line direction.
    slip_hkl : array-like object
        The Miller plane indices associated with the dislocation's slip plane.
    m : array-like object, optional
        The m unit vector for the solution.  m, n, and u (dislocation line)
        should be right-hand orthogonal.  Default value is [1, 0, 0]
        (x-axis).
    n : array-like object, optional
        The n unit vector for the solution.  m, n, and u (dislocation line)
        should be right-hand orthogonal.  Default value is [0, 1, 0]
        (y-axis).
    box : atomman.Box, optional
        The unit cell's box that the crystal vectors are taken with respect to. If
        not given, will use a cubic box with a=1 (vects are taken as Cartesian).
    tol : float, optional
        Tolerance parameter used to round off near-zero values.  Default
        value is 1e-8.
    
    Returns
    -------
    numpy.ndarray
        The (3,3) transformation matrix.
    """
    # Set default box to be cubic
    if box is None:
        box = Box()

    # Check m value
    if m is None:
        m = np.array([1,0,0], dtype=float)
    m = np.asarray(m, dtype=float)
    assert m.shape == (3, ), "m must be a 3D vector"
    assert np.isclose(np.linalg.norm(m), 1.0, atol = tol), "m must be a unit vector"

    # Check n value
    if n is None:
        n = np.array([0,1,0], dtype=float)
    n = np.asarray(n, dtype=float)
    assert n.shape == (3, ), "n must be a 3D vector"
    assert np.isclose(np.linalg.norm(n), 1.0, atol = tol), "n must be a unit vector"
    assert np.isclose(np.dot(m, n), 0.0, atol = tol), "m and n must be perpendicular"
    
    # ξ_axis is the Cartesian unit vector along dislocation line 
    ξ_axis = miller.vector_crystal_to_cartesian(ξ_uvw, box)
    ξ_axis = ξ_axis / np.linalg.norm(ξ_axis)
    
    # n_axis is the unit vector normal to slip plane
    n_axis = miller.plane_crystal_to_cartesian(slip_hkl, box)
    
    # m_axis is the unit vector in the "edge" direction
    m_axis = np.cross(n_axis, ξ_axis)
    
    transform = np.array([m_axis, n_axis, ξ_axis])
    
    # Transform transform to m,n,ξ orientation
    T = np.array([m,n,np.cross(m, n)]).T
    
    transform = T.dot(transform)
    
    return transform