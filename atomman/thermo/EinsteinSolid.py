# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
import atomman.unitconvert as uc

class EinsteinSolid():
    """
    Provides thermodynamic model predictions for an Einstein solid.
    """
    def __init__(self, theta: float, H0: float):
        """
        Class initializer

        Parameters
        ----------
        theta : float
            Einstein temperature in K.
        H0 : float
            The 0K enthalpy energy of the structure in eV/atom.  Identical
            to the cohesive energy of the structure at P = 0, T = 0K.
        """
        self.theta = theta
        self.H0 = H0

    @property
    def theta(self) -> float:
        """float: The Einstein temperature in K."""
        return self.__theta

    @theta.setter
    def theta(self, value: float):
        value = float(value)
        assert value > 0.0
        self.__theta = value

    @property
    def H0(self) -> float:
        """float: The 0K enthalpy energy in eV/atom."""
        return self.__H0

    @H0.setter
    def H0(self, value: float):
        self.__H0 = float(value)

    def H(self, T: npt.ArrayLike) -> np.ndarray:
        """
        The enthalpy as a function of T for the Einstein solid.

        H = H0 + (3 kB θ) / (exp(θ / T) - 1)
        
        Parameters
        ----------
        T : array-like object
            The temperatures at which to evaluate the enthalpy.

        Returns
        -------
        H : numpy.ndarray
            The enthalpy values at each given temperature.
        """

        θ = self.theta
        H0 = self.H0
        kB = uc.unit['kB']

        def ein(T):
            """Einstein model for H"""
            return H0 + (3 * kB * θ) / (np.exp(θ / T) - 1)
        
        def zero(T):
            """Return H0 for T=0 as above model is undef"""
            return H0
        
        return np.piecewise(T, [T > 0, T <= 0], [ein, zero])

    def Cv(self, T):
        """
        The volumetric heat capacity as a function of T for the Einstein solid.
        
        cV = (3 kB (θ / T)^2 exp(θ / T)) / (exp(θ / T) - 1)^2

        Parameters
        ----------
        T : array-like object
            The temperatures at which to evaluate the enthalpy.

        Returns
        -------
        Cv : numpy.ndarray
            The volumetric heat capacity values at each given temperature.
        """
        
        θ = self.theta
        kB = uc.unit['kB']
        
        def ein(T):
            """Einstein model for cV"""
            x = θ / T
            return 3 * kB * x**2 * np.exp(x) / (np.exp(x) - 1)**2
        
        def zero(T):
            """Return 0.0 for T=0 as above model is undef"""
            return 0

        return np.piecewise(T, [T > 0, T <= 0], [ein, zero])

    def G(self, T):
        """
        The Gibbs free energy as a function of T for the Einstein solid.
        
        G = H0 + 3/2 kB θ + 3 kB T ln(1 - exp(- θ / T))

        Parameters
        ----------
        T : array-like object
            The temperatures at which to evaluate the enthalpy.

        Returns
        -------
        G : numpy.ndarray
            The Gibbs free energy values at each given temperature.
        """
        θ = self.theta
        H0 = self.H0
        kB = uc.unit['kB']
        
        def ein(T):
            """Einstein model for G"""
            return H0 + 1.5 * kB * θ + 3 * kB * T * np.log(1 - np.exp(- θ / T))
        def zero(T):
            """Return H0 + 3/2 kB θ for T=0 as above model is undef"""
            return H0 + 1.5 * kB * θ
        
        return np.piecewise(T, [T > 0, T <= 0], [ein, zero])