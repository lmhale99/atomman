from copy import deepcopy

from scipy.spatial import Voronoi, ConvexHull

import numpy as np

from .. import Atoms

class InterstitialSite():
    
    def __init__(self,
                 pos: np.ndarray,
                 neighbor_atoms: Atoms):
        """
        Initializes an InterstitialSite object
        
        Parameters
        ----------
        pos : numpy.ndarray
            The coordinates for the interstitial position
        neighbor_pos : numpy.ndarray
            The coordinates for the atoms that neighbor the interstitial.
            These should correspond to either direct or replica atoms from
            the source atomic system.        
        """
        self.__pos = pos
        self.__neighbor_atoms = neighbor_atoms
    
    def __eq__(self, other):
        """
        Compare InterstitialSites using is_similar() with the default settings.
        """
        return self.is_similar(other)
        
        
    def is_similar(self,
                   other,
                   decimals: int = 6,
                   use_dmag: bool = False) -> bool:
        """
        Compares agains another InterstitialSite based on the neighbor atoms'
        atype and dvect or dmag values.
        
        
        Parameters
        ----------
        other : InterstitialSite
            The other InterstitialSite to compare against.
        decimals : int, optional
            The number of decimal points to round the dvect or dmag values to
            before comparing.  Default value is 6.
        use_dmag : bool
            If True then the comparison will use atype and dmag.  If False
            (default) then the comparison will use atype and dvect.  Roughly,
            this means that setting this to True will perform a rotation
            invariant comparison while leaving it False will not.
        """
        # First check number of atoms
        if self.neighbor_atoms.natoms != other.neighbor_atoms.natoms:
            return False
        
        atype0 = self.neighbor_atoms.atype
        atype1 = other.neighbor_atoms.atype

        if use_dmag:
            # Extract and round dmag values
            d0 = np.round(self.neighbor_dmag, decimals=decimals)
            d1 = np.round(other.neighbor_dmag, decimals=decimals)
        
            # Get sorting indices
            sorted_indices0 = np.lexsort((d0, atype0))
            sorted_indices1 = np.lexsort((d1, atype1))
            
        else:
            # Extract and round dvect values
            d0 = np.round(self.neighbor_dvect, decimals=decimals)
            d1 = np.round(other.neighbor_dvect, decimals=decimals)
        
            # Get sorting indices
            sorted_indices0 = np.lexsort((d0[:, 2], d0[:, 1], d0[:, 0], atype0))
            sorted_indices1 = np.lexsort((d1[:, 2], d1[:, 1], d1[:, 0], atype1))
        
        # Sort the arrays
        atype0 = atype0[sorted_indices0]
        atype1 = atype1[sorted_indices1]
        d0 = d0[sorted_indices0]
        d1 = d1[sorted_indices1]
            
        # Compare
        return np.allclose(atype0, atype1) and np.allclose(d0, d1)
        
    @property
    def pos(self):
        return self.__pos
    
    @property
    def neighbor_atoms(self):
        return self.__neighbor_atoms
    
    @property
    def volume(self):
        return ConvexHull(self.neighbor_atoms.pos).volume
    
    @property
    def neighbor_dvect(self):
        return self.neighbor_atoms.pos - self.pos
    
    @property
    def neighbor_dmag(self):
        return np.linalg.norm(self.neighbor_dvect, axis=1)
    
    def is_strained(self, rtol=1e-05, atol=1e-08):
        return not np.allclose(self.neighbor_dmag, self.neighbor_dmag[0])
    
    def asdict(self, rtol=1e-05, atol=1e-08):
        d = {
            'pos[0]': self.pos[0],
            'pos[1]': self.pos[1],
            'pos[2]': self.pos[2],
            '#neighbors': self.neighbor_atoms.natoms,
            'volume': self.volume,
            'strained': self.is_strained(rtol=rtol, atol=atol)
        }
        return d
    
def interstitial_site_finder(system):
    """
    Generates a list of interstitial sites for an atomic configuration using
    a Voronoi analysis.
    
    Parameters
    ----------
    system : atomman.System
        The atomic configuration to search for interstitial sites.
    
    Returns
    -------
    list of atomman.defect.InterstitialSite
        The identified interstitial sites.
    """
        
    # Supersize the system in all directions 
    bigsystem = system.supersize((-1,2), (-1,2), (-1,2))
    
    # Compute the Voronoi analysis
    vor = Voronoi(bigsystem.atoms.pos)
    
    # Filter out the Voronoi vertices that are not in the middle replica
    isin_ucell = system.box.inside(vor.vertices + .0005)
    vertices_pos = vor.vertices[isin_ucell]
    vertices_ids = np.where(isin_ucell)[0] # ids of vor.vertices that correspond to vertices
    
    # Initialize neighbors lists
    neighborlists = [ [] for nix in range(np.sum(isin_ucell))]
    
    # Search for all atoms that neighbor the vertices
    for atom_id, region_id in enumerate(vor.point_region):
        region_vertices_ids = vor.regions[region_id]
        for vertex_index, vertex_id in enumerate(vertices_ids):
            if vertex_id in region_vertices_ids:
                neighborlists[vertex_index].append(atom_id)
                
    interstitialsites = []
    for vertex_pos, neighborlist in zip(vertices_pos, neighborlists):
        neighbor_atoms = bigsystem.atoms[neighborlist]
        interstitialsites.append(InterstitialSite(vertex_pos, neighbor_atoms))
        
    return interstitialsites