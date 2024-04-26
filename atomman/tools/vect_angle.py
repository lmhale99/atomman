# coding: utf-8
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

def vect_angle(vect1: npt.ArrayLike,
               vect2: npt.ArrayLike,
               unit: str = 'degree') -> Union[float, np.ndarray]:
    """
    Returns the angle(s) between two sets of vectors.
    
    Parameters
    ----------
    vect1 : array-like object
        First vector or set of vectors.
    vect2 : array-like object
        Second vector or set of vectors.
    unit : str
        Specifies unit of returned angle: 'degree' or 'radian'. Default value
        is 'degree'.
        
    Returns
    -------
    float or numpy.NDArray
        The angle(s) between vect1 and vect2.
    """
    # Convert to numpy arrays if needed
    vect1 = np.asarray(vect1)
    vect2 = np.asarray(vect2)
    
    # Compute unit vectors
    u_vect1 = (vect1.T / np.linalg.norm(vect1, axis=-1)).T
    u_vect2 = (vect2.T / np.linalg.norm(vect2, axis=-1)).T
    
    # calculate cosine angle between vectors
    #cosine = np.dot(u_vect1, u_vect2)
    cosine = np.einsum('...i,...i', u_vect1, u_vect2)
    
    # Clean up any invalid numbers due to rounding
    try:
        # For multiple values
        cosine[cosine < -1] = -1
        cosine[cosine > 1] = 1
    except TypeError:
        # For single values
        if cosine < -1:
            cosine = -1
        elif cosine > 1:
            cosine = 1
    
    # Convert to degrees and return
    if unit == 'degree':
        return 180 * np.arccos(cosine) / np.pi
    elif unit == 'radian':
        return np.arccos(cosine)
    else:
        raise ValueError("unit must be 'degree' or 'radian'.")