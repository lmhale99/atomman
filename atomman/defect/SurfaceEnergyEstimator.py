# coding: utf-8
from typing import Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from ..tools import miller

class SurfaceEnergyEstimator():
    """
    Class that provides estimates of free surface energies for any miller index
    (hkl) of an FCC or BCC structure using the model provided in Seoane, A., &
    Bai, X. M. "A New Analytical Surface Energy Model for Arbitrary (hkl) Planes
    in BCC and FCC Metals". Surfaces and Interfaces, 103841 (2024).  The model
    uses the three basic surface energies of the (100), (110) and (111) planes
    as input.
    """
    def __init__(self,
                 e100: float,
                 e110: float,
                 e111: float,
                 structure: str):
        """
        Initialize the estimator by providing the three surface plane energy
        values and specifying the structure type.

        Parameters
        ----------
        e100 : float
            The surface energy of the (100) plane.
        e110 : float
            The surface energy of the (110) plane.
        e111 : float
            The surface energy of the (111) plane.
        structure : str
            Defines the crystal structure of the system.  This has to be either
            'fcc' or 'bcc'.
        """
        self.e100 = e100
        self.e110 = e110
        self.e111 = e111
        self.structure = structure

    @property
    def e100(self) -> float:
        """float: The surface energy of the (100) plane."""
        return self.__e100
    
    @e100.setter
    def e100(self, val: float):
        self.__e100 = float(val)

    @property
    def e110(self) -> float:
        """float: The surface energy of the (110) plane."""
        return self.__e110
    
    @e110.setter
    def e110(self, val: float):
        self.__e110 = float(val)

    @property
    def e111(self) -> float:
        """float: The surface energy of the (111) plane."""
        return self.__e111
    
    @e111.setter
    def e111(self, val: float):
        self.__e111 = float(val)

    @property
    def structure(self) -> str:
        """str: The crystal structure: 'fcc' or 'bcc'"""
        return self.__structure
    
    @structure.setter
    def structure(self, val: str):
        if val in ['fcc', 'bcc']:
            self.__structure = val
        else:
            raise ValueError("invalid structure value: must be 'fcc' or 'bcc'")
        
    def estimate(self,
                 hkl: Union[npt.ArrayLike, str]):
        """
        Estimate the surface energy for an arbitrary miller index (hkl).
        
        Parameters
        ----------
        hkl : array-like object or str
            The free surface plane to calculate expressed as integer
            Miller (hkl) indices.
        
        Returns
        -------
        surface_energy: float
            The surface energy of the desired (hkl) plane.
        """
        if self.structure == 'fcc':
            return self.estimate_fcc(hkl)
        elif self.structure == 'bcc':
            return self.estimate_bcc(hkl)
        else:
            raise ValueError('structure somehow not fcc or bcc!!!')

    def check_hkl(self,
                  hkl: Union[npt.ArrayLike, str]) -> np.ndarray:
        """
        Checks that a given (hkl) value is a 3-index integer Miller vector
        and if needed sorts it such that h >= k >= l.

        Parameters
        ----------
        hkl : array-like object or str
            The free surface plane to calculate expressed as integer
            Miller (hkl) indices.
        
        Returns
        -------
        hkl: numpy.ndarray
            The (hkl) plane with indices sorted h >= k >= l.
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

        return hkl
    
    def estimate_fcc(self, hkl):
        """
        Use the FCC structure formula from the paper

        α = |h - 2K + l| + |-h + 2K + l| + |-h + K + 2l|
        E_hkl = ( E_100 3(h - 3K - 2l + α) +
                  E_110 3sqrt(2)(3h + 3k - 2l - α) +
                  E_111 sqrt(3)(-3h + k + 10l + α) ) / (12 |hkl|)
        
        Parameters
        ----------
        hkl : array-like object or str
            The free surface plane to calculate expressed as integer
            Miller (hkl) indices.
        
        Returns
        -------
        surface_energy: float
            The surface energy of the desired (hkl) plane.
        """
        # Shortcut to the set energies
        e100 = self.e100
        e110 = self.e110
        e111 = self.e111
        
        # Check hkl value
        hkl = self.check_hkl(hkl)

        # Get the length of hkl
        length_hkl = np.sqrt(hkl.dot(hkl))

        # Define alpha FCC from the formula in the paper
        aFCC = (  abs(hkl[0] - 2 * hkl[1] + hkl[2]) 
                + abs(-hkl[0] + 2 * hkl[1] + hkl[2]) 
                + abs(-hkl[0] + hkl[1] + 2 * hkl[2]) )

        return (  3 * e100 * (hkl[0] - 3 * hkl[1] - 2 * hkl[2] + aFCC)
                + 3 * np.sqrt(2) * e110 * (3 * hkl[0] + 3 * hkl[1] - 2 * hkl[2] - aFCC)
                + np.sqrt(3) * e111 * (-3 * hkl[0] + hkl[1] + 10 * hkl[2] + aFCC)
                ) / (12 * length_hkl)
    
    def estimate_bcc(self, hkl):
        """
        Use the BCC structure formula from the paper

        α = |-h + K + l|
        E_hkl = ( E_100 3(h - K - l + α) +
                  E_110 3sqrt(2)(h + k - l - α) +
                  E_111 sqrt(3)(-h + k + 5l + α) ) / (6 |hkl|)
        
        Parameters
        ----------
        hkl : array-like object or str
            The free surface plane to calculate expressed as integer
            Miller (hkl) indices.
        
        Returns
        -------
        surface_energy: float
            The surface energy of the desired (hkl) plane.
        """
        # Shortcut to the set energies
        e100 = self.e100
        e110 = self.e110
        e111 = self.e111
        
        # Check hkl value
        hkl = self.check_hkl(hkl)

        # Get the length of hkl
        length_hkl = np.sqrt(hkl.dot(hkl))

        # Define alpha BCC from the formula in the paper
        aBCC = abs(-hkl[0] + hkl[1] + hkl[2])

        return (  3 * e100 * (hkl[0] - hkl[1] - hkl[2] + aBCC) 
                + 3 * np.sqrt(2) * e110 * (hkl[0] + hkl[1] - hkl[2] - aBCC)
                + np.sqrt(3) * e111 * (-hkl[0] + hkl[1] + 5 * hkl[2] + aBCC) 
               ) / (6 * length_hkl)