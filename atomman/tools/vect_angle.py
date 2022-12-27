# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

def vect_angle(vect1: npt.ArrayLike,
               vect2: npt.ArrayLike,
               unit: str = 'degree') -> float:
    """
    Returns the angle between two vectors.
    
    Parameters
    ----------
    vect1 : array-like object
        First vector.
    vect2 : array-like object
        Second vector.
    unit : str
        Specifies unit of returned angle: 'degree' or 'radian'. Default value
        is 'degree'.
        
    Returns
    -------
    float
        The angle between vect1 and vect2.
    """
    # Convert to numpy arrays if needed
    vect1 = np.asarray(vect1)
    vect2 = np.asarray(vect2)
    
    # Compute unit vectors
    u_vect1 = (vect1.T / np.linalg.norm(vect1, axis=-1)).T
    u_vect2 = (vect2.T / np.linalg.norm(vect2, axis=-1)).T
    
    # calculate cosine angle between vectors
    cosine = np.dot(u_vect1, u_vect2)
    
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