# coding: utf-8
# Standard Python libraries
import io
from typing import Optional, Tuple, Union
import uuid
from math import ceil

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_value, load_record
from yabadaba import unitconvert as uc

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from scipy.interpolate import CubicSpline

from matplotlib.figure import Figure
import matplotlib.pyplot as plt

__all__ = ['StructureFactor']

class StructureFactor(Record):

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'structure_factor'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'structure-factor'

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl html transformer"""
        return ('dbliquid.xsl', 'structure-factor.xsl')

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('dbliquid.xsd', 'structure-factor.xsd')


    ############################# Define Values  ##############################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """

        self._add_value('str', 'key', defaultvalue=str(uuid.uuid4()))
        self._add_value('str', 'id')
        self._add_value('longstr', 'source')
        self._add_value('str', 'element')
        self._add_value('float', 'temperature',
                        metadatakey='T (K)', unit='K')
        self._add_value('float', 'density', modelpath='atomic-density',
                        metadatakey='density (angstrom^-3)', unit='angstrom^-3')
        self._add_value('floatarray', 'q', valuerequired=True,
                        modelpath='structure-factor-plot.q', unit='angstrom^-1')
        self._add_value('floatarray', 'S', valuerequired=True,
                        modelpath='structure-factor-plot.S(q)')

    ########################## Additional settings ############################

    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return self.key

    ############################ Additional Methods ###########################

    def spline_S(self,
                 npoints: Optional[int] = None,
                 include_zero: bool = False,
                 inplace: bool = False):
        """
        Perform a cubic spline interpolation fit to the S(q) function.  One use
        for this is to interpolate missing intermediate points for S(q) tables
        that have irregularly spaced q values.

        Parameters
        ----------
        npoints : int or None, optional
            The number of q interpolation points to use in the range of q
            values.  If None (default) then the q values will be selected with
            a spacing equal to the smallest current q spacing.
        include_zero : bool, optional
            If True, then S(q) for q values less than given will be
            extrapolated by assuming S(q) smoothly goes to 0 at q=0.
        inplace: bool, optional
            If False (default), a new StructureFactor object will be returned
            with the splined values.  If True, the values of the current
            object will be updated instead.

        Returns
        -------
        StructureFactor
            If inplace is False, a new StructureFactor object containing the
            interpolated values will be returned.
        """

        if include_zero:
            # Add S(q=0) = 0
            q = np.array([0.0] + self.q.tolist())
            S = np.array([0.0] + self.S.tolist())
            bc_type = ['clamped', 'not-a-knot']

        else:
            q = self.q
            S = self.S
            bc_type = ['not-a-knot', 'not-a-knot']

        # Fit so S =0 and dS/dq = 0 at q = 0
        spline = CubicSpline(q, S, bc_type=bc_type)

        if npoints is not None:
            if include_zero:
                newq = np.linspace(0.0, q[-1], npoints+1)[1:]
            else:
                newq = np.linspace(q[0], q[-1], npoints)

        else:
            Δq = np.min(q[1:] - q[:-1])
            if include_zero:
                newq = np.arange(Δq, q[-1] + Δq, Δq)
            else:
                newq = np.arange(q[0], q[-1] + Δq, Δq)

        newS = spline(newq)

        if inplace:
            self.q = newq
            self.S = newS

        else:
            return StructureFactor(q=newq, S=newS, density=self.density,
                                   element=self.element, temperature=self.temperature,
                                   source = self.source)

    def rdf(self, zero_low_r: bool = False) -> Record:
        """
        Numerically transform S(q) into the radial distribution function, g(r).
        Note that there might be some numerical artifacts at low r values due
        to the maximum q being finite.
                
        Parameters
        ----------
        zero_low_r : bool
            The conversion can lead to some artifact fluctuations in g(r) at
            low r.  Setting this to True will cause all g(r) values to be zero
            below the largest r where g(r) is negative.  Default value is
            False.

        Returns
        -------
        RDF
            An RDF Record object with the estimated g(r) values.
        
        Raises
        ------
        ValueError
            If the q values are not evenly spaced.  Try using spline_S_q first.
        """
        if self.density is None:
            raise ValueError('density has not been set!')
        
        q = uc.get_in_units(self.q, 'Å^-1')
        S = self.S
        ρ = uc.get_in_units(self.density, 'Å^-3')
        π = np.pi

        Δq = q[1:] - q[:-1]
        if np.allclose(Δq, Δq[0]):
            Δq = Δq[0]
        else:
            raise ValueError('q must be evenly spaced to convert')

        # Compute the integral sums
        sums = np.sum(q * (S - 1) * np.sin(np.outer(q, q)), axis=1)

        # Convert q to r and compute g
        r = uc.set_in_units(q, 'Å')
        g = 1 + 1 / (2 * π**2 * ρ) * sums * Δq / q

        # Remove low r values once g = 0.0
        if zero_low_r:
            try:
                zeroi = np.where(g < 0.0)[0][-1] + 1
            except IndexError:
                pass
            else:
                g[:zeroi] = 0.0

        rdf = load_record('rdf', r=r, g=g, density=self.density,
                          element=self.element, temperature=self.temperature,
                          source = self.source)

        return rdf
    
    @staticmethod
    def hard_sphere_normalized(qr, packing):
        """
        Computes the structure factor for a single element hard-sphere fluid using the
        Percus-Yevick approximation and normalized coordinates.
        
        Parameters
        ----------
        qr : array-like
            Normalized qr coordinates: q multiplied by the hard-sphere radius, r
        packing : float
            The packing density. 
            
        Returns
        -------
        S: numpy.ndarray
            The hard-sphere structure factor for all q values.
        """
        
        # Shorthand eta to η
        η = packing
        
        # Percus-Yevick hard-sphere parameters
        a = (1 + 2 * η)**2 / (1 - η)**4
        b = -6 * η * (1 + η / 2)**2 / (1 - η)**4
        c = η / 2 * (1 + 2 * η)**2 / (1 - η)**4

        # x is qd
        x = 2 * qr
        g = (a / x**2 * (np.sin(x) - x * np.cos(x))
        + b / x**3 * (2 * x * np.sin(x) + (2 - x**2) * np.cos(x) - 2)
        + c / x**5 * (-x**4 * np.cos(x) + 4 * ( (3 * x**2 - 6) * np.cos(x) + x * (x**2 - 6) * np.sin(x) + 6) ) )
            
        return 1 / (1 + 24 * η * g / x)

    @classmethod
    def hard_sphere(cls, q, radius, density):
        """
        Computes the structure factor for a single element hard-sphere fluid using the
        Percus-Yevick approximation.
        
        Parameters
        ----------
        q : array-like
            The q coordinates to compute the structure factor at.
        radius : float
            The hard-sphere radius.
        density : float
            The atomic number density.
            
        Returns
        -------
        StructureFactor
            The computed hard sphere structure factor.
        """
        # Define equation terms
        π = np.pi
        d = 2 * radius
        ρ = density
        
        # Compute the packing factor
        η = π * ρ * d**3 / 6
        
        qr = q * radius
        
        S = cls.hard_sphere_normalized(qr, η)

        return cls(q=q, S=S, density=density)
    
    def S_q_plot(self,
                length_unit: str = 'Å',
                figsize: Optional[tuple] = None,
                fig: Optional[Figure] = None,
                **kwargs) -> Figure:
        """
        Convenience method for generating g(r) plots.

        Parameters
        ----------
        length_unit : str, optional
            The unit of length to use with the q values: q values will be
            displayed in 1/length_units.  Default value is Å for angstrom.
        figsize : tuple, optional
            The x,y size of the figure to return.  Default value is (10, 6).
        fig : matplotlib.figure, optional
            An existing figure object to add the new plot to.  If not given, a
            new figure is generated.
        **kwargs : dict, optional
            Additional keywords are passed into the underlying 
            matplotlib.pyplot.plot() call. This allows control of such things
            like line color, style, etc.
        """
        
        # Get values and convert units
        inv_length_unit = f'{length_unit}^-1'
        q = uc.get_in_units(self.q, inv_length_unit)
        S = self.S
        qmax = ceil(q.max())
        Smax = ceil(S.max())

        # Create new figure
        if fig is None:
            if figsize is None:
                figsize = (10, 6)
            fig = plt.figure(figsize=figsize)

        # Check bounds of existing figure
        else:
            old_qmax = fig.axes[0].get_xlim()[-1]
            old_Smax = fig.axes[0].get_ylim()[-1]
            if old_qmax > qmax:
                qmax = old_qmax
            if old_Smax > Smax:
                Smax = old_Smax
        
        # Plot with/without formatting options
        if 'fmt' in kwargs:
            fmt = kwargs.pop('fmt')
            plt.plot(q, S, fmt, **kwargs)
        else:
            plt.plot(q, S, **kwargs)
        
        # Set labels and ranges
        plt.xlabel(f'q ({inv_length_unit})', size='x-large')
        plt.ylabel('S(q)', size='x-large')
        plt.xlim(0.0, qmax)
        plt.ylim(0.0, Smax)

        return fig