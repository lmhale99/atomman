# coding: utf-8

# Standard Python libraries
from typing import Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
import atomman.unitconvert as uc
from ..tools import aslist

class IdealGas():
    
    def __init__(self,
                 T: float,
                 V: float,
                 m: Union[float, npt.ArrayLike],
                 N: Union[float, npt.ArrayLike]):
        """
        Parameters
        ----------
        T : float
            The temperature to use.
        V : float
            The total volume of the system.
        m : float or list
            The atomic masses for each atom type.
        N : float or list
            The number of atoms for each atom type.
        """
        
        self.__T = T
        self.__V = V
        self.__m = np.array(aslist(m))
        self.__N = np.array(aslist(N))
        
    @property
    def T(self):
        return self.__T
    
    @property
    def V(self):
        return self.__V
    
    @property
    def m(self):
        return self.__m
    
    @property
    def N(self):
        return self.__N

    @property
    def rho(self):
        return self.N.sum() / self.V
    
    @property
    def Λ(self):

        T = self.T
        m = self.m
        π = np.pi
        h = uc.unit['hPlanck']
        kB = uc.unit['kB']
        
        # Compute de Broglie thermal wavelength(s)
        return (h**2 / (2 * π * kB * T * m)) ** 0.5
    
    @property
    def c(self):
        return self.N / self.N.sum()
    
    @property
    def S(self):
        
        Λ = self.Λ
        ρ = self.rho
        c = self.c
        kB = uc.unit['kB']
        
        return self.N.sum() * kB * (5 / 2 - (c * np.log(ρ * c * Λ**3)).sum() )
    
    @property
    def F(self):
        
        T = self.T
        N = self.N
        Λ = self.Λ
        ρ = self.rho
        c = self.c
        kB = uc.unit['kB']

        return (N * kB * T * (np.log(ρ) + 3 * np.log(Λ) - 1 + np.log(c))).sum()