#Standard library imports
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
from . import Atoms
from . import Box
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
            self.atoms_prop('pos', atoms.pos, scale=True)
    
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
    
    def atoms_prop(self, key, value=None, scale=False):
        """
        Wrapper around atoms.prop dictionary adding a scale argument.
        
        Parameters
        ----------
        key : str
            The property key name to get/set values for.
        value : numpy.ndarray, optional
            The value to assign to key.
        scale : bool, optional
            Flag indicating if values should be scaled/unscaled. If True:
            - Values being retrieved are scaled from absolute Cartesian to box relative vectors.
            - Values being set are unscaled from box relative to absolute Cartesian vectors.
            Default value is False.
        """
        
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Get value
        if value is None:
            value = self.atoms.prop[key]
            if scale is True:
                value = self.scale(value)
            return value
        
        # Set value
        else:
            if scale is True:
                value = self.unscale(value)
            self.atoms.prop[key] = value
    
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
            self.atoms_prop('pos', spos, scale=True)
        
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
        self.atoms_prop('pos', spos, scale=True)
        
        # Modify box vectors and origin by new min and max
        origin = self.box.origin + mins.dot(self.box.vects) 
        avect = self.box.avect * (maxs[0] - mins[0])
        bvect = self.box.bvect * (maxs[1] - mins[1])
        cvect = self.box.cvect * (maxs[2] - mins[2])
        self.box_set(avect=avect, bvect=bvect, cvect=cvect, origin=origin)