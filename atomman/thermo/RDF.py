# coding: utf-8

# https://pandas.pydata.org/
import pandas as pd

# http://www.numpy.org/
import numpy as np

# Local imports
import atomman.unitconvert as uc
from ..tools import uber_open_rmode

class RDF():
    """Class for interpreting LAMMPS rdf (radial distribution function) output files"""
    
    def __init__(self, rdf_file=None):
        
        if rdf_file is not None:
            self.read(rdf_file)
            
    def read(self, rdf_file): 
        
        self.data = {}
        numvals = 0
        timeheader = 3
        with uber_open_rmode(rdf_file) as rdf_data:
            
            for i, line in enumerate(rdf_data):
                line = line.decode('UTF-8')
                    
                if i == 0:
                    self.info = line[2:]
                elif i < 3:
                    continue
                
                elif i == timeheader:
                    terms = line.split()
                    timestep = int(terms[0])
                    numvals = int(terms[1])
                    
                    # initialize data for the timestep
                    self.data[timestep] = pd.DataFrame({'r': np.empty(numvals),
                                                        'g': np.empty(numvals),
                                                        'coord': np.empty(numvals)})
                    timeheader += numvals + 1
                    
                else:
                    terms = line.split()
                    j = int(terms[0]) - 1
                    self.data[timestep].loc[j, 'r'] = float(terms[1])
                    self.data[timestep].loc[j, 'g'] = float(terms[2])
                    self.data[timestep].loc[j, 'coord'] = float(terms[3])
    
    @property
    def r(self):
        """
        The r coordinates for the first measurement - assuming all measurements are the same
        """
        return self.data[list(self.data.keys())[0]].r


    @property
    def g(self):
        """
        The average g values at all r across all measurements
        """
        g = np.zeros_like(self.r)
        count = 0
        for timestep in self.data:
            g += self.data[timestep].g
            count += 1
        
        return g / count

    @property
    def I(self):
        """
        The integrand used to estimate 2-body excess entropy from g(r)
        """
        def ln(x):
            x = np.asarray(x, dtype=float)
            lnx = np.log(x)
            lnx[np.isneginf(lnx)] = 0.0
            return lnx

        r = self.r
        g = self.g
        
        return (g * ln(g) - g + 1) * r**2

    def S2body(self, density: float) -> float:
        """
        The 2-body excess entropy estimated from g(r)

        Parameters
        ----------
        density : float
            The atomic density

        Returns
        -------
        float
            The 2-body excess entropy per atom.
        """
        π = np.pi
        kB = uc.unit['kB']
        ρ = density

        return -2 * π * ρ * kB * np.trapz(self.I, self.r)