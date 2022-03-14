# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

def axes_check(axes: npt.ArrayLike, 
               tol: float = 1e-8) -> np.ndarray:
    """
    Checks that given axes are orthogonal and right-handed.
    
    Parameters
    ----------
    axes : array-like object
        A 3x3 array of axes vectors.
    tol : float, optional
        Tolerance to use in checking if axes are orthogonal and right-handed.
        Default value is 1e-8.
    
    Returns
    -------
    numpy.ndarray
        A 3x3 array of the corresponding unit vectors of axes.
    
    Raises
    ------
    ValueError
        If axes are not orthogonal or right-handed.
    """
    # Convert to numpy.ndarray if needed
    axes = np.asarray(axes)
    assert axes.shape == (3,3), 'Invalid axes shape'
    
    # Compute unit axes
    uaxes = (axes.T / np.linalg.norm(axes, axis=1)).T
    
    # Check if orthogonal
    if not np.allclose(np.dot(uaxes, uaxes.T), np.identity(3), atol=tol):
        raise ValueError('axes are not orthogonal')
    
    # Check if right-handed
    if not np.allclose(np.cross(uaxes[0], uaxes[1]), uaxes[2], atol=tol):
        raise ValueError('axes are not right-handed')
    
    return uaxes