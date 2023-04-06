# coding: utf-8
import io
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
import atomman.unitconvert as uc
from ..tools import uber_open_rmode

class RDF():
    """Handles tabular radial distribution function data and computes derived properties"""
    
    def __init__(self,
                 r: Optional[npt.ArrayLike] = None,
                 g: Optional[npt.ArrayLike] = None,
                 coord: Optional[npt.ArrayLike] = None,
                 lammps_rdf_file: Union[str, io.IOBase, None] = None,
                 density: Optional[float] = None):
        """
        Class init.
        
        Parameters
        ----------
        r: array-like object, optional
            r coordinates.
        g: array-like object, optional
            g(r) values.
        coord: array-like object, optional
            coord(r) values.
        rdf_file : str, file-like object, optional
            A LAMMPS rdf output file to extract values of r, g, and coord.
            Cannot be given with r, g, or coord parameters.
        density : float, optional
            The atomic number density of the system.  This is required to
            calculate some derived properties.
        """
        # load from lammps rdf file if given
        if lammps_rdf_file is not None:
            if r is not None or g is not None or coord is not None:
                raise ValueError('lammps_rdf_file cannot be given with r, g, or coord')
            self.read_lammps_file(lammps_rdf_file)
        
        # Set values if directly given
        else:
            self.r = r
            self.g = g
            self.coord = coord
        
        # Set density
        self.density = density

        # Set derived properties
        self.__I = None
        self.__structure_factor = None
        self.__entropy_2body = None

    def clear_derived_properties(self):
        """Resets values of cached properties back to None"""
        self.__I = None
        self.__structure_factor = None
        self.__entropy_2body = None

    @property
    def r(self) -> np.ndarray:
        """numpy.ndarray: The r coordinates"""
        if self.__r is not None:
            return self.__r
        else:
            raise ValueError('r has not been set')

    @r.setter
    def r(self, value: npt.ArrayLike):
        if value is None:
            self.__r = None
        else:
            self.__r = np.asarray(value, dtype=float)
        self.clear_derived_properties()

    @property
    def g(self) -> np.ndarray:
        """numpy.ndarray: The g(r) values"""
        if self.__g is not None:
            return self.__g
        else:
            raise ValueError('g has not been set')

    @g.setter
    def g(self, value: npt.ArrayLike):
        if value is None:
            self.__g = None
        else:
            self.__g = np.asarray(value, dtype=float)
        self.clear_derived_properties()

    @property
    def coord(self) -> np.ndarray:
        """numpy.ndarray: The coord(r) values"""
        if self.__coord is not None:
            return self.__coord
        else:
            raise ValueError('coord has not been set')

    @coord.setter
    def coord(self, value: npt.ArrayLike):
        if value is None:
            self.__coord = None
        else:
            self.__coord = np.asarray(value, dtype=float)
        self.clear_derived_properties()

    @property
    def density(self) -> float:
        """float: The atomic number density"""
        if self.__density is not None:
            return self.__density
        else:
            raise ValueError('density has not been set')
        
    @density.setter
    def density(self, value: float):
        if value is None:
            self.__density = None
        else:
            self.__density = float(value)
        self.clear_derived_properties()
    
    @property
    def I(self) -> np.ndarray:
        """numpy.ndarray: The integrand used to estimate 2-body excess entropy from g(r)"""
        if self.__I is None:
        
            def ln(x):
                """natural log modification to handle invalid numbers"""
                return np.piecewise(x, [x > 0], [np.log, lambda x: 0])

            r = self.r
            g = self.g
            
            self.__I = (g * ln(g) - g + 1) * r**2

        return self.__I

    @property
    def structure_factor(self) -> np.ndarray:
        """numpy.ndarray: The structure factor calculated from g(r)"""
        if self.__structure_factor is None:

            r = self.r
            g = self.g
            ρ = self.density
            π = np.pi
            
            dr = r[1] - r[0]
        
            # Compute the integral sums
            sums = np.sum(r * (g - 1) * np.sin(np.outer(r, r)), axis=1)
        
            self.__structure_factor = 1 + 4 * π * ρ * sums * dr / r

        return self.__structure_factor

    @property
    def entropy_2body(self) -> float:
        """float: The 2-body excess entropy estimated from g(r)"""
        if self.__entropy_2body is None:
            π = np.pi
            kB = uc.unit['kB']
            ρ = self.density

            self.__entropy_2body = -2 * π * ρ * kB * np.trapz(self.I, self.r)
        
        return self.__entropy_2body

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
