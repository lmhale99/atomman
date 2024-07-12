# coding: utf-8
from typing import Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from ..tools import miller

def surface_energy_estimate(hkl: Union[npt.ArrayLike, str],
                            e100: float,
                            e110: float,
                            e111: float,
                            structure: str) -> float:
    """
    Calculates the surface energy of any miller index (hkl) for any FCC
    or BCC structure after using the model provided in Seoane, A., & Bai, 
    X. M. "A New Analytical Surface Energy Model for Arbitrary (hkl) 
    Planes in BCC and FCC Metals". Surfaces and Interfaces, 103841 (2024).
    The model uses the three basic surface energies of the (100), (110) and 
    (111) planes as input.
    
    Parameters
    ----------
    hkl : array-like object or str
        The free surface plane to to calculate expressed as integer
        Miller (hkl) indices.
    e100 : float
        The surface energy of the (100) plane.
    e110 : float
        The surface energy of the (110) plane.
    e111 : float
        The surface energy of the (111) plane.
    structure : str
        It defines if the system is 'fcc' or 'bcc'.
        This has to be either 'fcc' or 'bcc'.
      
    Returns
    -------
    surface_energy : float
        The surface energy of the desired hkl.
    """
    # Read string hkl
    if isinstance(hkl, str):
        hkl = miller.fromstring(hkl)

    # Check hkl values
    hkl = np.asarray(hkl)
    if hkl.shape != (3,):
        raise ValueError('invalid hkl indices: must be 3 values')
    
    if np.allclose(hkl, np.asarray(hkl, dtype=int)):
        hkl = np.asarray(hkl, dtype=int)
    else:
        raise ValueError('hkl indices must be integers')
    
    if np.allclose(hkl, np.zeros(3)):
        raise ValueError('hkl indices cannot all be zero')

    # Rearrange so h >= k >= l
    hkl = np.flip(np.sort(hkl))
    
    # Get the length of hkl
    length_hkl = np.sqrt(hkl.dot(hkl))
    
    # Check the structure and apply the corresponding formula
    if structure == 'bcc':
        # Define alpha BCC from the formula in the paper
        aBCC = abs(-hkl[0] + hkl[1] + hkl[2])

        return (  3 * e100 * (hkl[0] - hkl[1] - hkl[2] + aBCC) 
                + 3 * np.sqrt(2) * e110 * (hkl[0] + hkl[1] - hkl[2] - aBCC)
                + np.sqrt(3) * e111 * (-hkl[0] + hkl[1] + 5 * hkl[2] + aBCC) 
               ) / (6 * length_hkl)

    elif structure == 'fcc':
        # Define alpha FCC from the formula in the paper
        aFCC = (  abs(hkl[0] - 2 * hkl[1] + hkl[2]) 
                + abs(-hkl[0] + 2 * hkl[1] + hkl[2]) 
                + abs(-hkl[0] + hkl[1] + 2 * hkl[2]) )

        return (  3 * e100 * (hkl[0] - 3 * hkl[1] - 2 * hkl[2] + aFCC)
                + 3 * np.sqrt(2) * e110 * (3 * hkl[0] + 3 * hkl[1] - 2 * hkl[2] - aFCC)
                + np.sqrt(3) * e111 * (-3 * hkl[0] + hkl[1] + 10 * hkl[2] + aFCC)
                ) / (12 * length_hkl)
    
    else:
        raise ValueError('The structure has not been defined as either fcc or bcc')
