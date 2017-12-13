# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from . import Atoms
from . import Box
import atommantest.core.supersize
import atommantest.core.dvect
import atommantest.core.NeighborList
import atommantest.core.load
import atommantest.convert.system_model
from ..compatibility import iteritems, range

class System(object):
    """
    A representation of an atomic system.  This combines Atoms with Box and
    adds methods and attributes involving both.
    """
    
    def __init__(self, atoms=Atoms(), box=Box(), pbc=(True, True, True), scale=False):
        """
        Initilize a System by joining an am.Atoms and am.Box instance.
        
        Parameters
        ----------
        atoms : atomman.Atoms
            The underlying Atoms object to build system around.
        box : atomman.Box
            The underlying box object to build system around.
        pbc : tuple or list of bool
            Indicates which of the dimensions related to the three box vectors
            are periodic.  Default value is (True, True, True).
        scale : bool
            If True, atoms.pos will be scaled relative to the box.  Default
            value is False.
        """
        
        # Check data types
        if not isinstance(atoms, Atoms):
            raise TypeError('Invalid atoms type')
        if not isinstance(box, Box):
            raise TypeError('Invalid box type')
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Set properties
        self.__atoms = atoms
        self.__box = box
        self.pbc = pbc

        # Rescale pos if needed
        if scale is True:
            self.atoms_prop('pos', value=atoms.pos, scale=True)
    
    def __str__(self):
        """str : The string representation of a system."""
        return '\n'.join([str(self.box),
                          'natoms = ' + str(self.natoms),
                          'natypes = ' + str(self.natypes),
                          str(self.atoms)])
    
    def __len__(self):
        return self.natoms
    
    @property
    def atoms(self):
        """atomman.Atoms : underlying Atoms object."""
        return self.__atoms
    
    @property
    def natoms(self):
        """int : The number of atoms in the Atoms class."""
        return self.__atoms.natoms
    
    @property
    def atypes(self):
        """list of int : List of unique atom types."""
        return self.__atoms.atypes
    
    @property
    def natypes(self):
        """int : The number of atom types in the Atoms class."""
        return self.__atoms.natypes
    
    @property
    def box(self):
        """atommman.Box : underlying Box object."""
        return self.__box
    
    @property
    def pbc(self):
        """list of bool : The periodic boundary condition settings."""
        return self.__pbc
        
    @pbc.setter
    def pbc(self, value):
        pbc = np.asarray(value, dtype=bool)
        assert pbc.shape == (3,), 'invalid pbc entry' 
        self.__pbc = pbc
    
    def atoms_prop(self, key=None, index=None, value=None, a_id=None, scale=False):
        """
        Extends Atoms.prop() by adding a scale argument.
        
        Parameters
        ----------
        key : str, optional
            Per-atom property name.
        index : int, list, slice, optional
            Index of atoms.
        value : any, optional
            Property values to assign.
        a_id : int, optional
            Integer atom index.  Left in for backwards compatibility.
        scale : bool, optional
            Flag indicating if values should be scaled/unscaled. If True:
            - Values being retrieved are scaled from absolute Cartesian to box relative vectors.
            - Values being set are unscaled from box relative to absolute Cartesian vectors.
            Default value is False.
        """
        
        # Check that scale is bool
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Call atoms.prop() for scale is False
        if scale is False:
            if value is None:
                return self.atoms.prop(key=key, index=index, a_id=a_id)
            else:
                self.atoms.prop(key=key, index=index, value=value, a_id=a_id)
        
        # Handle property scaling if scale is True
        else:
            # Handle a_id
            if a_id is not None:
                if index is not None:
                    raise ValueError('a_id and index cannot both be given')
                index = a_id
            
            # Get values if value is None
            if value is None:
                # If no key, return copy of atoms (all or a slice) with scaled pos
                if key is None:
                    if index is None:
                        newatoms = deepcopy(self.atoms)
                    else:
                        newatoms = deepcopy(self.atoms[index])
                    newatoms.pos = self.scale(newatoms.pos)
                    return newatoms
                
                # If key is given, return scaled property values
                else:
                    if index is None:
                        value = self.atoms.view[key]
                    else:
                        value = self.atoms.view[key][index]
                    return self.scale(value)
            
            # Set values if value is given
            else:
                # If no key, unscale pos of Atoms value and set to atoms (or a slice)
                if key is None:
                    if not isinstance(value, Atoms):
                        raise TypeError('If key is None, value must be instance of atomman.Atoms')
                    
                    value.pos = self.unscale(value.pos)
                    if index is None:
                        self.atoms[:] = value
                    else:
                        self.atoms[index] = value
                
                # If key is given, unscale value and set to property
                else:
                    value = self.unscale(value)
                    if index is None:
                        self.atoms.view[key] = value
                    else:
                        self.atoms.view[key][index] = value
    
    def box_set(self, **kwargs):
        """
        Extends box.set() with a scale argument.
        
        Parameters
        ----------
        scale : bool
            If True, the scaled (box-relative) positions remain unchanged.
            Default value is False (absolute positions unchanged).
        other kwargs : any
            Any other kwargs allowed by box.set().
        """
        
        # Pop scale
        scale = kwargs.pop('scale', False)
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Hold scaled positions constant
        if scale is True:
            spos = self.atoms_prop('pos', scale=True)
            self.box.set(**kwargs)
            self.atoms_prop('pos', value=spos, scale=True)
        
        # Call box.set without scaling
        else:
            self.box.set(**kwargs)
    
    def scale(self, value):
        """
        Scales 3D vectors from absolute Cartesian coordinates to relative box
        coordinates.
        
        Parameters
        ----------
        value : numpy.ndarray
            Values to scale.
        """
        
        # Retrieve parameters
        value = np.asarray(value)
        vects = self.box.vects
        inverse = np.linalg.inv(vects)
        origin = self.box.origin
        
        # Convert
        return (value - origin).dot(inverse)
          
    def unscale(self, value):
        """
        Unscales 3D vectors from relative box coordinates to absolute
        Cartesian coordinates.
        
        Parameters
        ----------
        value : numpy.ndarray
            Values to unscale.
        """
        
        # Retrieve parameters
        value = np.asarray(value)
        vects = self.box.vects
        origin = self.box.origin
        
        # Convert
        return value.dot(vects) + origin
        
    def wrap(self):
        """
        Wrap atoms around periodic boundaries and extend non-periodic
        boundaries such that all atoms are within the box.
        """
        
        # mins and maxs are box dimensions relative to box vectors, i.e 0 to 1
        mins = np.array([0.0, 0.0, 0.0])
        maxs = np.array([1.0, 1.0, 1.0])
        
        # Retrieve scaled pos
        spos = self.atoms_prop('pos', scale=True)
        
        # Loop over three pbc directions
        for i in range(3):
            
            # Wrap atoms across periodic boundaries
            if self.pbc[i]:
                spos[:, i] -= np.floor(spos[:, i])
            
            # Shift min and max to encompass atoms across non-periodic bounds
            else:
                min = spos[:, i].min()
                max = spos[:, i].max()
                if min < mins[i]: mins[i] = min - 0.001
                if max > maxs[i]: maxs[i] = max + 0.001
        
        # Unscale spos and save to pos
        self.atoms_prop('pos', value=spos, scale=True)
        
        # Modify box vectors and origin by new min and max
        origin = self.box.origin + mins.dot(self.box.vects) 
        avect = self.box.avect * (maxs[0] - mins[0])
        bvect = self.box.bvect * (maxs[1] - mins[1])
        cvect = self.box.cvect * (maxs[2] - mins[2])
        self.box_set(avect=avect, bvect=bvect, cvect=cvect, origin=origin)
    
    def dvect(self, pos_0, pos_1, code=None):
        """
        Computes the shortest distance between pos_0 and pos_1 using box 
        dimensions and accounting for periodic boundaries.
        
        Parameters
        ----------
        pos_0 : numpy.ndarray or index
            Absolute Cartesian vector position(s) to use as reference point(s).
            If the value can be used as an index, then self.atoms.pos[pos_0]
            is taken.
        pos_1 : numpy.ndarray or index
            Absolute Cartesian vector position(s) to find relative to pos_0.
            If the value can be used as an index, then self.atoms.pos[pos_1]
            is taken.
        code: str, optional
            Option for specifying which underlying code function to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        """
        # Test if pos_0 and pos_1 can be used as numpy array indices
        try:
            pos_0 = self.atoms.view['pos'][pos_0]
        except:
            pos_0 = np.asarray(pos_0)
        try:
            pos_1 = self.atoms.view['pos'][pos_1]
        except:
            pos_1 = np.asarray(pos_1)
        
        # Call atomman.core.dvect using self's box and pbc
        return atommantest.core.dvect(pos_0, pos_1, self.box, self.pbc, code=code)
    
    def neighborlist(self, **kwargs):
        """
        Builds a neighbor list for the system.  The resulting NeighborList
        object is saved to the object as attribute 'neighbors'.
        
        Parameters
        ----------
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Must be given if
            model is not given.
        model : str or file-like object, optional
            Gives the file path or content to load.  If given, initialsize is
            the only other allowed parameter.
        cmult : int, optional
            Parameter associated with the binning routine.  Default value is most
            likely the fastest.
        code: str, optional
            Option for specifying which underlying code function of nlist to use:
            - 'cython' uses the version of the function built in cython (faster).
            - 'python' uses the purely python version.
            Default is 'cython' if the code can be imported, otherwise 'python'.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        """
        if 'system' in kwargs:
            raise KeywordError("Cannot give 'system' as it is taken as the current object")
        else:
            kwargs['system'] = self
        self.neighbors = atommantest.core.NeighborList(**kwargs)
    
    def supersize(self, a_size, b_size, c_size):
        """
        Increases the system's size by replicating it along the box vectors.
        
        The multiplier values *_size are taken to be integer tuples (m, n) where
        m <= 0 and n >= 0.  The system multiplication works such that if n = -m,
        then the seed system's origin will be at the center of the new system.
        If only one integer is given, then it is assigned to m or n depending on
        its sign, and the other value is taken to be 0.
        
        Parameters
        ----------
        a_size : int or tuple of int
            Single int or two integers specifying replication along the avect
            direction.
        b_size -- int or tuple of int
            Single int or two integers specifying replication along the bvect
            direction.
        c_size -- int or tuple of int
            Single int or two integers specifying replication along the cvect
            direction.
        """
        system = atommantest.core.supersize(self, a_size, b_size, c_size)
        self.__box = system.box
        self.__atoms = system.atoms
        self.pbc = system.pbc
    
    def load(self, style, input, **kwargs):
        """
        Read in system content from another format.  The system will
        be updated and the list of atomic symbols added as a symbols attribute.
        
        Parameters
        ----------
        style : str
            Indicates the format of the content to load as an atomman.System
        input : str, file-like object or object
            The content to load.
        kwargs
            Any extra keyword arguments to pass to the underlying load methods.
        """
        system, symbols = atommantest.core.load(style, input, **kwargs)
        self.__box = system.box
        self.__atoms = system.atoms
        self.pbc = system.pbc
        self.symbols = symbols
    
    def model(self, **kwargs):
        """
        Return a DataModelDict 'cell' representation of the system.
        
        Parameters
        ----------
        box_unit : str, optional
            Length unit to use for the box. Default value is 'Angstrom'.
        symbols : list, optional
            list of atom-model symbols corresponding to the atom types.
        elements : list, optional
            list of element tags corresponding to the atom types.
        prop_units : dict, optional
            dictionary where the keys are the property keys to include, and
            the values are units to use. If not given, only the positions in
            scaled units are included.
        a_std : float, optional
            Standard deviation of a lattice constant to include if available.
        b_std : float, optional
            Standard deviation of b lattice constant to include if available.
        c_std : float, optional
            Standard deviation of c lattice constant to include if available.
        """
            
        return atommantest.convert.system_model.dump(self, **kwargs)