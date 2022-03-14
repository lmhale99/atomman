# coding: utf-8

# Standard Python imports
from __future__ import annotations
from typing import Callable, Optional, Union

#https://matplotlib.org/
import matplotlib.pyplot as plt

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from . import gradient, integrator
import atomman.unitconvert as uc

class BasePath():
    """
    Generic class representing an energy path.
    """
    
    def __init__(self,
                 coord: npt.ArrayLike,
                 energyfxn: Callable,
                 gradientfxn: str = 'cdiff',
                 gradientkwargs: Optional[dict] = None,
                 integratorfxn: str = 'rk'):
        """
        Class initializer.
        
        Parameters
        ----------
        coord : array-like object
            The list of coordinates associated with the points along the path.
        energyfxn : function
            The function that evaluates the energy associated with the different
            point coordinates.
        gradientfxn : str or function, optional
            The function to use to estimate the gradient of the energy.  A str
            value of 'central_difference' or 'cdiff' (default) will use
            atomman.mep.gradient.central_difference.
        gradientkwargs : dict, optional
            The keyword arguments (i.e. settings) to use with the gradientfxn.
            Default is an empty dictionary, i.e. default settings of gradientfxn.
        integratorfxn : str or function, optional
            The function to use to integrate relaxation steps.  A str value of
            'euler' will use atomman.mep.integrator.euler, while a str value of
            'rungekutta' or 'rk' (default) will use atomman.mep.integrator.rungekutta.
        """
        
        if isinstance(coord, BasePath):
            coord = coord.coord
        self.coord = coord

        if callable(energyfxn):
            self.__energyfxn = energyfxn
        else:
            raise TypeError('energyfxn must be a callable object')
        
        self.gradientfxn = gradientfxn
        self.integratorfxn = integratorfxn
        
        if gradientkwargs is None:
            self.__gradientkwargs = {}
        if isinstance(gradientkwargs, dict):
            self.__gradientkwargs = gradientkwargs
        else:
            raise TypeError('gradientkwargs must be None or a dict')
    
    @property
    def coord(self) -> np.ndarray:
        """numpy.NDArray : The coordinates for each point along the path."""
        return self.__coord
    
    @coord.setter
    def coord(self, value: npt.ArrayLike):
        self.__coord = np.asarray(value)
    
    @property
    def arccoord(self) -> np.ndarray:
        """numpy.NDArray : The arc length coordinates for each point along the path."""
        
        # Compute the arc length distance between each point
        s = np.zeros(len(self.coord))
        s[1:] = np.linalg.norm(self.coord[1:] - self.coord[:-1], axis=1)
        
        # Compute the arc length coordinates
        α = np.empty(len(self.coord))
        for i in range(len(α)):
            α[i] = s[:i+1].sum()
        
        return α

    @property
    def energyfxn(self) -> Callable:
        """function : The function for evaluating the energy"""
        return self.__energyfxn
    
    @property
    def gradientfxn(self) -> Callable:
        """function : The function for evaluating the energy gradients"""
        return self.__gradientfxn
    
    @gradientfxn.setter
    def gradientfxn(self, value: Union[str, Callable]):
        if isinstance(value, str):
            if value == 'central_difference' or value == 'cdiff':
                self.__gradientfxn = gradient.central_difference
            else:
                raise ValueError('Unknown gradientfxn style')
        elif callable(value):
            self.__gradientfxn = value
        else:
            raise TypeError('gradientfxn must be a str of a known style or a callable object')

    @property
    def gradientkwargs(self) -> dict:
        """dict : The keyword arguments to use when calling gradientfxn"""
        return self.__gradientkwargs
    
    @property
    def integratorfxn(self) -> Callable:
        """function : The function to use for integrating minimization steps"""
        return self.__integratorfxn
    
    @integratorfxn.setter
    def integratorfxn(self, value: Union[str, Callable]):
        if isinstance(value, str):
            if value == 'rungekutta' or value == 'rk':
                self.__integratorfxn = integrator.rungekutta
            elif value == 'euler':
                self.__integratorfxn = integrator.euler
            else:
                raise ValueError('Unknown integratorfxn style')
        elif callable(value):
            self.__integratorfxn = value
        else:
            raise TypeError('integratorfxn must be a str of a known style or a callable object')

    @property
    def unittangent(self) -> np.ndarray:
        """numpy.NDArray : The tangent vectors along the path at each point."""
        raise NotImplementedError('Defined by subclasses')
    
    def energy(self,
               coord: Optional[npt.ArrayLike] = None) -> np.ndarray:
        """
        Evaluates energy values associated with either the path points or
        given coordinates.
        
        Parameters
        ----------
        coord : array-like object, optional
            Coordinates to evaluate the energies at.  If not given, will use
            the object's coord values.
            
        Returns
        -------
        np.NDArray
            The evaluated energies.
        """
        if coord is None:
             coord = self.coord
        
        return self.energyfxn(coord)
    
    def grad_energy(self,
                    coord: Optional[npt.ArrayLike] = None) -> np.ndarray:
        """
        Evaluates the gradient of the energy associated with either the path
        points or given coordinates.
        
        Parameters
        ----------
        coord : array-like object, optional
            Coordinates to evaluate the energy gradients at.  If not given,
            will use the object's coord values.
            
        Returns
        -------
        np.NDArray
            The evaluated energies.
        """
        if coord is None:
             coord = self.coord
        return self.gradientfxn(self.energyfxn, coord, **self.gradientkwargs)
    
    @property
    def force(self) -> np.ndarray:
        """numpy.NDArray : The computed force associated with moving along the path at each point."""
        return np.einsum('ij,ij->i', self.grad_energy(), self.unittangent)

    def step(self, *args, **kwargs) -> BasePath:
        """
        Performs a single relaxation step.
        
        Returns
        -------
        newpath : BasePath
            A Path with coordinates evolved forward by one timestep.
        """
        raise NotImplementedError('Defined by subclasses')

    def relax(self, *args, **kwargs) -> BasePath:
        """
        Perform multiple relaxation and/or climb steps until either the
        maximum coordinate displacement per step drops below a tolerance or
        the maximum number of steps is reached.

        Returns
        -------
        newpath : BasePath
            A Path with coordinates evolved forward from the relaxation.
        """
        raise NotImplementedError('Defined by subclasses')

    def plot_energy(self,
                    energy_unit: Optional[str] = None,
                    length_unit: Optional[str] = None,
                    ax: Optional[plt.axes] = None,
                    **kwargs) -> Optional[plt.figure]:
        """
        Creates a plot of the energies along the path as a function of the
        arc coordinates.

        Parameters
        ----------
        energy_unit : str or None, optional
            If given, the energy values will be converted from atomman's
            working units to the specified units.  Default value of None will
            do no conversions, which is useful if the energyfxn is not reporting
            in atomman's working units.
        length_unit : str or None, optional
            If given, the arc coordinates will be converted from atomman's
            working units to the specified units.  Default value of None will
            do no conversions, which is useful if the coords are not given
            in atomman's working units.
        ax : matplotlib.pyplot.axis
            A pre-existing plotting axis.  Allows for more control and the
            use of subplots.
        **kwargs : any, optional
            All additional keyword arguments will be passed to 
            matplotlib.pyplot.figure().

        Returns
        -------
        matplotlib.pyplot.figure
            The generated figure allowing for further modifications.  Returned
            if ax is None.
        """
        if ax is None:
            fig = plt.figure(**kwargs)
            ax = fig.add_subplot(111)
            returnfig = True
        else:
            assert len(kwargs) == 0, 'ax and extra kwargs cannot both be given'
            returnfig = False

        # Handle energy units
        energy = self.energy()
        if energy_unit is None:
            ax.set_ylabel('Energy', size='x-large')
        else:
            ax.set_ylabel(f'Energy ({energy_unit})', size='x-large')
            energy = uc.get_in_units(energy, energy_unit)

        # Handle length unit
        arccoord = self.arccoord
        if length_unit is None:
            ax.set_xlabel('Arc coordinate', size='x-large')
        else:
            ax.set_xlabel(f'Arc coordinate ({length_unit})', size='x-large')
            arccoord = uc.get_in_units(arccoord, length_unit)
        ax.set_xlim(0, arccoord[-1])

        ax.plot(arccoord, energy, lw=3)

        if returnfig:
            return fig

    def plot_force(self,
                   force_unit: Optional[str] = None,
                   length_unit: Optional[str] = None,
                   ax: Optional[plt.axes] = None,
                   **kwargs) -> Optional[plt.figure]:
        """
        Creates a plot of the force to move along the path as a function of the
        arc coordinates.

        Parameters
        ----------
        force_unit : str or None, optional
            If given, the force values will be converted from atomman's
            working units to the specified units.  Default value of None will
            do no conversions, which is useful if the energyfxn is not reporting
            in atomman's working units.
        length_unit : str or None, optional
            If given, the arc coordinates will be converted from atomman's
            working units to the specified units.  Default value of None will
            do no conversions, which is useful if the coords are not given
            in atomman's working units.
        ax : matplotlib.pyplot.axis
            A pre-existing plotting axis.  Allows for more control and the
            use of subplots.
        **kwargs : any, optional
            All additional keyword arguments will be passed to 
            matplotlib.pyplot.figure().

        Returns
        -------
        matplotlib.pyplot.figure
            The generated figure allowing for further modifications.  Returned
            if ax is None.
        """
        if ax is None:
            fig = plt.figure(**kwargs)
            ax = fig.add_subplot(111)
            returnfig = True
        else:
            assert len(kwargs) == 0, 'ax and extra kwargs cannot both be given'
            returnfig = False

        # Handle force units
        force = self.force
        if force_unit is None:
            ax.set_ylabel('Force', size='x-large')
        else:
            ax.set_ylabel(f'Force ({force_unit})', size='x-large')
            force = uc.get_in_units(force, force_unit)

        # Handle length unit
        arccoord = self.arccoord
        if length_unit is None:
            ax.set_xlabel('Arc coordinate', size='x-large')
        else:
            ax.set_xlabel(f'Arc coordinate ({length_unit})', size='x-large')
            arccoord = uc.get_in_units(arccoord, length_unit)
        ax.set_xlim(0, arccoord[-1])

        ax.plot(arccoord, force, lw=3)

        if returnfig:
            return fig