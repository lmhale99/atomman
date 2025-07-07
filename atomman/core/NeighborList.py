# coding: utf-8

# Standard Python imports
import io
from typing import Union, Optional

# http://www.numpy.org/
import numpy as np

# atomman imports
from .nlist import nlist
from ..tools import uber_open_rmode

class NeighborList(object):
    """Class that finds and stores the neighbor atoms for a system."""
    
    def __init__(self, **kwargs):
        """
        Class initializer.  Calls NeighborList.load() if
        'model' is given, otherwise calls NeighborList.build().
        
        Parameters
        ----------
        system : atomman.System, optional
            The system to calculate the neighbor list for. Must be given if
            model is not given.
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Must be given if
            model is not given.
        model : str or file-like object, optional
            Gives the file path or content to load.  If given, no other
            parameters are allowed.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
        nlist : array-like object, optional
            An array consisting of int values, where each row is associated with
            an atom, the first value in each row is that atom's coordination, and
            all subsequent values (up to the coordination number) are indices of
            the neighboring atoms.
        """
        # Load from model
        if 'model' in kwargs:
            model = kwargs.pop('model')
            self.load(model, **kwargs)
        
        # Convert from existing nlist array
        elif 'nlist' in kwargs:
            nlist = kwargs.pop('nlist')
            assert len(kwargs) == 0
            self.__nlist = nlist
            self.__coord = nlist[:, 0]
            self.__neighbors = nlist[:, 1:]
        
        # Build new neighbor list
        else:
            system = kwargs.pop('system')
            cutoff = kwargs.pop('cutoff')
            self.build(system, cutoff, **kwargs)


    @classmethod
    def freud(cls,
              system,
              cutoff: Optional[float] = None,
              coord: Optional[int] = None,
              return_bonds: bool = False):
        """
        Use the AABBQuery from the freud Python package to build the neighbor
        list.  This has better performance than the native atomman method and
        allows for searching by number of neighbors, but requires additional
        package installs. 

        Parameters
        ----------
        system : atomman.System
            The system to calculate the neighbor list for. 
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Either cutoff
            or coord are required.
        coord : int, optional
            Number of neighbors to find for each atom.  Either cutoff
            or coord are required.
        return_bonds : bool, optional
            If True, then the method will also return the list of bonds
            found by the AABBQuery.  Useful if you also want the dmag vects
            for all neighbor pairs.
        """
        try:
            import freud
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError('package freud must be installed for this method') from e
        
        # Build freud query terms dict
        if cutoff is not None:
            if coord is not None:
                raise ValueError('cutoff and coord cannot both be given')
            query_terms = dict(
                mode = 'ball',
                r_max = cutoff,
                exclude_ii = True)
            
        elif coord is not None:
            query_terms = dict(
                mode = 'nearest',
                num_neighbors = coord,
                exclude_ii = True)
        
        else:
            raise ValueError('either cutoff or coord must be given')

        # Convert atomman system into freud box and points
        box, points = system.dump('freud')

        # Build query object and perform the query
        aq = freud.locality.AABBQuery(box, points)      
        bonds = aq.query(points, query_terms)
        coord = 15
        # Find max coord
        if coord is None:
            coord = 0       # coord is max coord across all atoms
            last_i = -1
            c = 0           # c is current coord count for atom i
            for bond in bonds:
                i, j, dmag = bond
                if i != last_i:
                    if c > coord:
                        coord = c
                    last_i = i
                    c = 1
                else:
                    c += 1

        # Build nlist
        nlist = np.empty((system.natoms, coord+1), dtype=np.int64)
        nlist[:, 0] = 0

        for bond in bonds:
            i, j, dmag = bond
            
            # Increase atom's coordination
            nlist[i, 0] += 1
            
            # Add j to the nlist row
            nlist[i, nlist[i, 0]] = j

        # Initialize neighborlist object
        neighbors = cls(nlist=nlist)

        if return_bonds:
            return neighbors, bonds
        else:
            return neighbors
        

    @property
    def coord(self) -> np.ndarray:
        """numpy.ndarray: The atomic coordination numbers"""
        return self.__coord
    
    @property
    def nlist(self) -> np.ndarray:
        """numpy.ndarray: The underlying numpy array of coord + neighbor ids"""
        return self.__nlist
    
    def __len__(self) -> int:
        """len returns the number of atoms"""
        return len(self.__coord)
    
    def __getitem__(self, key):
        """Get returns the list of neighbors for the specified atom."""
        return self.__neighbors[key, :self.coord[key]]
    
    def build(self,
              system,
              cutoff: float,
              initialsize: int = 20,
              deltasize: int = 10):
        """
        Builds the neighbor list for a system.
        
        Parameters
        ----------
        system : atomman.System 
            The system to calculate the neighbor list for.
        cutoff : float
            Radial cutoff distance for identifying neighbors.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
        """
        # Call nlist
        self.__nlist = nlist(system, cutoff, initialsize=initialsize, 
                             deltasize=deltasize)
        
        # Split coord and neighbors
        self.__coord = self.__nlist[:, 0]
        self.__neighbors = self.__nlist[:, 1:]
        
    def load(self, model: Union[str, io.IOBase]):
        """
        Read in a neighbor list from a file.
        
        Parameters
        ----------
        model : str or file-like object
            Gives the file path or content to load.
        """
        # First pass determines number of atoms and max number of neighbors
        nterms = 0
        natoms = 0
        with uber_open_rmode(model) as fin:
            for line in fin:
                line = line.decode('UTF-8')
                terms = line.split()
                n_n = len(terms)
                if terms[0][0] != '#' and n_n > 0:
                    natoms += 1
                    if n_n > nterms:
                        nterms = n_n
            
            self.__nlist = np.empty((natoms, nterms+1), dtype=int)
            self.__coord = self.__nlist[:, 0]
            self.__neighbors = self.__nlist[:, 1:]
            
            # Second pass gets values
            self.__coord[:] = 0
            fin.seek(0)
            for line in fin:
                line = line.decode('UTF-8')
                terms = line.split()
                if len(terms) > 0 and terms[0][0] != '#':
                    i = int(terms[0])
                    self.__coord[i] = len(terms) - 1
                    for j in range(1, len(terms)):
                        self.__neighbors[i, j-1] = terms[j]
    
    def dump(self, fname: str):
        """
        Saves the neighbor list to a file.
        
        Parameters
        ----------
        fname : str
            The file name to save the content to.
        """
        with open(fname, 'w') as fp:
            fp.write('# Neighbor list:\n')
            fp.write('# The first column gives an atom index.\n')
            fp.write('# The rest of the columns are the indexes of the identified neighbors.\n')
            for i in range(len(self)):
                fp.write('%i' % i)
                for j in self[i]:
                    fp.write(' %i' % j)
                fp.write('\n')