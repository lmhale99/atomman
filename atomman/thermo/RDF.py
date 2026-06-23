# Standard Python libraries
import io
from typing import Optional, Tuple, Union
import uuid
from math import ceil

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_record
from yabadaba import unitconvert as uc

# http://www.numpy.org/
import numpy as np

from scipy.interpolate import CubicSpline

from matplotlib.figure import Figure
import matplotlib.pyplot as plt


# Local imports
from ..tools import uber_open_rmode

__all__ = ['RDF']

class RDF(Record):
    """Handles tabular radial distribution function data and computes derived properties"""

    def __init__(self,
                 model: Union[str, io.IOBase, DM, None] = None,
                 name: Optional[str] = None,
                 **kwargs):
        """
        Initializes an RDF Record object.
        
        Parameters
        ----------
        model : str, file-like object or DataModelDict, optional
            A JSON/XML data model for the content.
        name : str, optional
            The unique name to assign to the record.  If model is a file
            path, then the default record name is the file name without
            extension.
        r: array-like object, optional
            The radial coordinates.
        g: array-like object, optional
            The radial distribution function, g(r), values.
        coord: array-like object, optional
            Coordination, coord(r), values.
        density : str, optional
            The atomic density.  This is required to calculate some derived
            properties.
        temperature : str, optional
            The temperature in Kelvin.
        element : str, optional
            The chemical element symbol.
        source : str, optional
            Information on the record's source, i.e. citation and/or creator
            details.
        lammps_rdf_file : str, file-like object, optional
            A LAMMPS rdf output file to extract values of r, g, and coord.
            Cannot be given with r, g, or coord parameters.
        """
        lammps_rdf_file = None
        if 'lammps_rdf_file' in kwargs:
            if model is not None:
                raise ValueError('model and lammps_rdf_file cannot both be given')
            if 'r' in kwargs or 'g' in kwargs or 'coord' in kwargs:
                raise ValueError('lammps_rdf_file cannot be given with r, g, or coord')
            lammps_rdf_file = kwargs.pop('lammps_rdf_file')

        super().__init__(model=model, name=name, **kwargs)
        
        # load from lammps rdf file if given
        if lammps_rdf_file is not None:
            self.read_lammps_file(lammps_rdf_file)

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'rdf'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'rdf'

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl html transformer"""
        return ('dbliquid.xsl', 'rdf.xsl')

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('dbliquid.xsd', 'rdf.xsd')
    
    ############################# Define Values  ##############################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """

        self._add_value('str', 'key', 
                                defaultvalue=str(uuid.uuid4()))
        self._add_value('str', 'id')
        self._add_value('longstr', 'source')
        self._add_value('str', 'relax_liquid_key')
        self._add_value('str', 'element')
        self._add_value('float', 'temperature',
                        metadatakey='T (K)', unit='K')
        self._add_value('float', 'density', modelpath='atomic-density',
                        metadatakey='density (angstrom^-3)', unit='angstrom^-3')
        self._add_value('floatarray', 'r', valuerequired=True,
                        modelpath='radial-density-function-plot.r',
                        unit='angstrom')
        self._add_value('floatarray', 'g', valuerequired=True,
                        modelpath='radial-density-function-plot.g(r)')
        self._add_value('floatarray', 'coord',
                        modelpath='radial-density-function-plot.coord(r)')
                                  
    ########################## Additional settings ############################
    
    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return self.key

    ############################ Additional Methods ###########################

    @classmethod
    def from_system(cls,
                    system,
                    nbins: int = 400,
                    rmin: float = 0.0,
                    rmax: float = 10.0,
                    **kwargs):
        """
        Compute the radial distribution function for an atomman system.

        Parameters
        ----------
        system : atomman.System
            The system to evaluate.
        nbins : int
            The number of bins to use.  Default is 400.
        rmin : float
            The minimum radial distance to include. Default is 0.0.
        rmax : float
            The maximum radial distance to include. Default is 10.0.
        **kwargs
            Any other keyword arguments recognized by the class initialization
            that do not conflict with the r, g, coord and density values
            computed here.
        
        Returns
        -------
        RDF : The radial distribution object.
        """
        # Get system information
        system.wrap()
        natoms = system.natoms
        volume = system.box.volume
        density = natoms / volume

        # Create supercell and nlist
        sizemults = [1, 1, 1]
        if system.pbc[0]:
            sizemults[0] = ceil(2 * rmax / system.box.lx)
        if system.pbc[1]:
            sizemults[1] = ceil(2 * rmax / system.box.ly)
        if system.pbc[2]:
            sizemults[2] = ceil(2 * rmax / system.box.lz)
        bigsystem = system.supersize(*sizemults)
        nlist = bigsystem.neighborlist(cutoff=rmax)

        # Create empty histogram and find the edges
        counts, edges = np.histogram(np.array([]), nbins, range=(rmin, rmax), density=False)
        mid_r = (edges[1:] + edges[:-1]) / 2
        delta_r = edges[1] - edges[0]

        # Compute the denominator
        denom = (4 * np.pi * natoms * (natoms - 1) / volume * mid_r**2) * delta_r

        # Loop over all atom ids in original system
        for i_id in range(natoms):
            
            # Find distances to the neighbor atoms
            dmags = bigsystem.dmag(system.atoms.pos[i_id], nlist[i_id])

            # Add distances to the histogram
            counts += np.histogram(dmags, nbins, range=(rmin, rmax), density=False)[0]

        # Remove 0 count
        if rmin == 0:
            counts[0] = 0
        
        # Compute coord and g
        coord = np.cumsum(counts) / natoms
        g = counts / denom

        # Build new RDF object
        return cls(r=mid_r, g=g, coord=coord, density=density, **kwargs)
        
    @classmethod
    def average(cls,
                rdfs: list,
                **kwargs):
        """
        Compiles an averaged RDF object from multiple other RDF objects.  The
        input RDFs should all have the same radial bins.  The final RDF will
        have average values for g, coord, and density.
        
        Parameters
        ----------
        rdfs : list
            A list of RDF objects to average over.
        **kwargs
            Any other keyword arguments recognized by the class initialization
            that do not conflict with the r, g, coord and density values
            computed here.
        """
        # Init values
        gs = []
        coords = []
        densities = []
        r = None
        
        # Loop over input rdfs
        for rdf in rdfs:
            
            # Check r values
            if r is None:
                r = rdf.r
            elif r.shape[0] != rdf.r.shape[0] or not np.allclose(r, rdf.r):
                raise ValueError('r values differ!')

            # Build list of terms
            gs.append(rdf.g)
            coords.append(rdf.coord)
            densities.append(rdf.density)

        # Average terms
        g = np.mean(gs, axis=0)
        coord = np.mean(coords, axis=0)
        density = np.mean(densities)

        # Build new RDF object
        return cls(r=r, g=g, coord=coord, density=density, **kwargs)

    def integrand(self) -> np.ndarray:
        """
        Convert g(r) to I(r), which is the integrand used to estimate 2-body
        excess entropy: I(r) = ( g(r) * ln(g(r)) - g(r) + 1 ) * r^2
        
        Returns
        -------
        numpy.ndarray
            The I(r) values
        """
        def ln(x):
            """natural log modification to handle invalid numbers"""
            return np.piecewise(x, [x > 0], [np.log, lambda x: 0])

        r = self.r
        g = self.g

        I = (g * ln(g) - g + 1) * r**2

        return I

    def spline_g(self,
                 npoints: Optional[int] = None,
                 include_zero: bool = False,
                 inplace: bool = False):
        """
        Perform a cubic spline interpolation fit to the g(r) function.  One use
        for this is to interpolate missing intermediate points for g(r) tables
        that have irregularly spaced r values.

        Parameters
        ----------
        npoints : int or None, optional
            The number of r interpolation points to use in the range of r
            values.  If None (default) then the r values will be selected with
            a spacing equal to the smallest current r spacing.
        include_zero : bool, optional
            If True, then g(r) for r values less than given will be
            set to 0.0.
        inplace: bool, optional
            If False (default), a new RDF object will be returned
            with the splined values.  If True, the values of the current
            object will be updated instead.

        Returns
        -------
        RDF
            If inplace is False, a new RDF object containing the
            interpolated values will be returned.
        """

        r = self.r
        g = self.g
        bc_type = ['not-a-knot', 'not-a-knot']
        spline = CubicSpline(r, g, bc_type=bc_type)

        if npoints is not None:
            if include_zero:
                newr = np.linspace(0.0, r[-1], npoints+1)[1:]
            else:
                newr = np.linspace(r[0], r[-1], npoints)

        else:
            Δr = np.min(r[1:] - r[:-1])
            if include_zero:
                newr = np.arange(Δr, r[-1] + Δr, Δr)
            else:
                newr = np.arange(r[0], r[-1] + Δr, Δr)

        newg = spline(newr)
        if include_zero:
            # Zero any low r g(r) values not given
            newg[newr < r[0]] = 0.0

            # Search for any g(r) < 0 and zero it and any g(r) for r < it.
            ltzero = np.where(newg < 0.0)[0]
            if len(ltzero) > 0:
                end = ltzero[-1] + 1
                newg[:end] = 0.0


        if inplace:
            self.r = newr
            self.g = newg

        else:
            return RDF(r=newr, g=newg, density=self.density,
                       element=self.element, temperature=self.temperature,
                       source = self.source)

    def structure_factor(self) -> Record:
        """
        Numerically transform g(r) into the structure factor, S(q).  Note that
        there might be some numerical artifacts at low q values due to the 
        maximum r being finite.

        Returns
        -------
        StructureFactor
            A StructureFactor Record object with the estimated S(q) values.
        """
        if self.density is None:
            raise ValueError('density has not been set!')

        r = uc.get_in_units(self.r, 'Å')
        g = self.g
        ρ = uc.get_in_units(self.density, 'Å^-3')
        π = np.pi

        Δr = r[1:] - r[:-1]
        if np.allclose(Δr, Δr[0]):
            Δr = Δr[0]
        else:
            raise ValueError('q must be evenly spaced to convert')

        # Compute the integral sums
        sums = np.sum(r * (g - 1) * np.sin(np.outer(r, r)), axis=1)

        # Convert r to q and compute S_q
        q = uc.set_in_units(r, 'Å^-1')
        S = 1 + 4 * π * ρ * sums * Δr / r

        sf = load_record('structure_factor', q=q, S=S, density=self.density,
                          element=self.element, temperature=self.temperature,
                          source = self.source)

        return sf

    def entropy_2body(self) -> float:
        """
        Estimate the 2-body excess entropy from g(r) by integrating I(r).
        
        Returns
        -------
        float
            The estimated 2-body excess entropy.
        """
        if self.density is None:
            raise ValueError('density has not been set!')
        
        π = np.pi
        kB = uc.unit['kB']
        ρ = self.density
        I = self.integrand()

        entropy_2body = -2 * π * ρ * kB * np.trapz(I, self.r)

        return entropy_2body

    @classmethod
    def from_relax_liquid(cls, record):
        """
        Extracts the RDF data from an iprPy relax_liquid calculation record and
        generates an rdf record.
        """
        f = record.get_file('rdf_dump.txt')
        rdf = cls(
            id = f'{record.system.composition}--{record.temperature:.2f}--md--{record.potential.potential_LAMMPS_id}',
            lammps_rdf_file = f,
            element = record.system.composition,
            temperature = record.temperature,
            relax_liquid_key = record.key,
            source = f'MD results from iprPy relax_liquid calculation {record.name}',
            density = 1 / record.volume)
        f.close()

        return rdf

    def read_lammps_file(self,
                         lammps_rdf_file: Union[str, io.IOBase]): 
        """
        Read a LAMMPS rdf output file and extract r, g(r) and coord(r) values.
        The r values are be taken from the first timestep, while the g(r) and
        coord(r) values are computed as the mean values at each r across all
        timesteps.
        
        Parameters
        ----------
        lammps_rdf_file : str, file-like object
            The LAMMPS rdf output file to read.
        """
        # Set default values
        numrows = 0
        headerindex = 3
        numtimesteps = 0
        rowindex = 0

        with uber_open_rmode(lammps_rdf_file) as rdf_data:

            for lineindex, line in enumerate(rdf_data):

                # Skip the first three lines
                if lineindex < 3:
                    continue

                line = line.decode('UTF-8')

                # Handle timestep headers
                if lineindex == headerindex:
                    terms = line.split()

                    # Initialize empty data on the first header
                    if lineindex == 3:
                        numrows = int(terms[1])
                        r = np.zeros(numrows, dtype=float)
                        g = np.zeros(numrows, dtype=float)
                        coord = np.zeros(numrows, dtype=float)
                        first = True

                    # Verify same number of rows for later series
                    elif numrows != int(terms[1]):
                        raise ValueError('g(r) tabulations not uniform for all timesteps!')

                    else:
                        first = False

                    # Reset rowindex
                    rowindex = 0

                    # Update headerindex and numtimesteps
                    headerindex += numrows + 1
                    numtimesteps += 1

                # Handle data lines
                else:
                    terms = line.split()

                    # Set r values only if first
                    if first:
                        r[rowindex] = float(terms[1])

                    # Add g, coord values
                    g[rowindex] += float(terms[2])
                    coord[rowindex] += float(terms[3])

                    # Update rowindex
                    rowindex += 1

        # Set class attributes
        self.r = r
        self.g = g / numtimesteps
        self.coord = coord / numtimesteps

    def build_lammps_file(self,
                          lammps_rdf_file: Optional[str] = None):

        nrows = self.r.shape[0]

        # Header lines
        lines = [
            '# RDF data computed from dump files',
            '# Unused Number-of-rows',
            '# Row r rdf coord',
            f'0 {nrows}'
        ]

        # Tabulation
        for i in range(nrows):
            lines.append(f'{i+1} {self.r[i]} {self.g[i]} {self.coord[i]}')

        contents = '\n'.join(lines)

        if lammps_rdf_file is None:
            return contents
        else:
            with open(lammps_rdf_file, 'w') as f:
                f.write(contents)

    def g_r_plot(self,
                length_unit: str = 'Å',
                figsize: Optional[tuple] = None,
                fig: Optional[Figure] = None,
                **kwargs) -> Figure:
        """
        Convenience method for generating g(r) plots.

        Parameters
        ----------
        length_unit : str, optional
            The unit of length to display the r values in.  Default value
            is Å for angstrom.
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
        r = uc.get_in_units(self.r, length_unit)
        g = self.g
        rmax = ceil(r.max())
        gmax = ceil(g.max())

        # Create new figure
        if fig is None:
            if figsize is None:
                figsize = (10, 6)
            fig = plt.figure(figsize=figsize)

        # Check bounds of existing figure
        else:
            old_rmax = fig.axes[0].get_xlim()[-1]
            old_gmax = fig.axes[0].get_ylim()[-1]
            if old_rmax > rmax:
                rmax = old_rmax
            if old_gmax > gmax:
                gmax = old_gmax
        
        # Plot with/without formatting options
        if 'fmt' in kwargs:
            fmt = kwargs.pop('fmt')
            plt.plot(r, g, fmt, **kwargs)
        else:
            plt.plot(r, g, **kwargs)
        
        # Set labels and ranges
        plt.xlabel(f'r ({length_unit})', size='x-large')
        plt.ylabel('g(r)', size='x-large')
        plt.xlim(0.0, rmax)
        plt.ylim(0.0, gmax)

        return fig
    
    def coord_r_plot(self,
                     length_unit: str = 'Å',
                     figsize: Optional[tuple] = None,
                     fig: Optional[Figure] = None,
                     **kwargs) -> Figure:
        """
        Convenience method for generating coord(r) plots.

        Parameters
        ----------
        length_unit : str, optional
            The unit of length to display the r values in.  Default value
            is Å for angstrom.
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
        r = uc.get_in_units(self.r, length_unit)
        coord = self.coord
        rmax = ceil(r.max())
        coordmax = ceil(coord.max())

        # Create new figure
        if fig is None:
            if figsize is None:
                figsize = (10, 6)
            fig = plt.figure(figsize=figsize)

        # Check bounds of existing figure
        else:
            old_rmax = fig.axes[0].get_xlim()[-1]
            old_coordmax = fig.axes[0].get_ylim()[-1]
            if old_rmax > rmax:
                rmax = old_rmax
            if old_coordmax > coordmax:
                coordmax = old_coordmax
        
        # Plot with/without formatting options
        if 'fmt' in kwargs:
            fmt = kwargs.pop('fmt')
            plt.plot(r, coord, fmt, **kwargs)
        else:
            plt.plot(r, coord, **kwargs)
        
        # Set labels and ranges
        plt.xlabel(f'r ({length_unit})', size='x-large')
        plt.ylabel('coord(r)', size='x-large')
        plt.xlim(0.0, rmax)
        plt.ylim(0.0, coordmax)

        return fig
    
    def integrand_r_plot(self,
                         length_unit: str = 'Å',
                         figsize: Optional[tuple] = None,
                         fig: Optional[Figure] = None,
                         **kwargs) -> Figure:
        """
        Convenience method for generating I(r) plots.

        Parameters
        ----------
        length_unit : str, optional
            The unit of length to display the r values in.  Default value
            is Å for angstrom.
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
        r = uc.get_in_units(self.r, length_unit)
        integrand = self.integrand()
        rmax = ceil(r.max())
        integrandmax = ceil(integrand.max())

        # Create new figure
        if fig is None:
            if figsize is None:
                figsize = (10, 6)
            fig = plt.figure(figsize=figsize)

        # Check bounds of existing figure
        else:
            old_rmax = fig.axes[0].get_xlim()[-1]
            old_integrandmax = fig.axes[0].get_ylim()[-1]
            if old_rmax > rmax:
                rmax = old_rmax
            if old_integrandmax > integrandmax:
                integrandmax = old_integrandmax
        
        # Plot with/without formatting options
        if 'fmt' in kwargs:
            fmt = kwargs.pop('fmt')
            plt.plot(r, integrand, fmt, **kwargs)
        else:
            plt.plot(r, integrand, **kwargs)
        
        # Set labels and ranges
        plt.xlabel(f'r ({length_unit})', size='x-large')
        plt.ylabel('I(r)', size='x-large')
        plt.xlim(0.0, rmax)
        plt.ylim(0.0, integrandmax)

        return fig