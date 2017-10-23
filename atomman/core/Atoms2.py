# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy
from collections import OrderedDict

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# atomman imports
from ..compatibility import iteritems, int, inttype, range

class Atoms(object):
    """
    Class for representing a collection of atoms.
    """
    
    class PropertyDict(OrderedDict):
        """Extends OrderedDict to work with Atoms"""
        
        def __init__(self, host):
            """PropertyDict needs to know its host Atoms object."""
            self.__host = host
            super(Atoms.PropertyDict, self).__init__()
        
        def __setitem__(self, key, value):
            """
            Modifies OrderedDict__setitem__() such that
            1. All values are converted to numpy.ndarrays.
            2. For new keys, checks are done that data is compatible with host.natoms.
            3. New keys are also added as attributes to host Atoms object, if allowed.
            4. For existing keys, new values are saved over the old ones as opposed to
               only changing the name assignment.

            Parameters
            ----------
            key : str
                The property key name to assign values to.
            value : numpy.ndarray
                The values to assign.
            """
            # Shortcut to host
            host = self.__host
            
            # Python 2: change to unicode if needed.
            try:
                key = key.decode('UTF-8')
            except:
                pass
            
            # Convert to numpy.ndarray if needed
            value = np.asarray(value)

            # Broadcast if needed and allowed
            if value.shape == ():
                value = np.array(np.broadcast_to(value, (host.natoms,) + value.shape))
            elif value.shape[0] == 1:
                value = np.array(np.broadcast_to(value, (host.natoms,) + value.shape[1:]))
            elif value.shape[0] != host.natoms:
                raise ValueError('First dimension of value must be 1 or natoms')

            # If key is already assigned, save value over existing values
            if key in self.keys():
                self[key][:] = value
            
            # Otherwise, set new item and try assigning attribute to host
            else:
                super(Atoms.PropertyDict, self).__setitem__(key, value)
                try:
                    assert key not in dir(host)
                    super(Atoms, host).__setattr__(key, value)
                except:
                    pass
    
    def __init__(self, natoms=None, atype=None, pos=None, **kwargs):
        """
        Class initializer.
        
        Parameters
        ==========
        natoms : int
            The number of atoms.
        atype : int or list/ndarray of int
            The integer atomic types to assign to all atoms.
        pos : list/ndarray of float
            The atomic positions to assign to all atoms.
        kwargs : dict
            All additional key/value pairs are assigned as properties.
            
        Returns
        =======
        Atoms
            The Atoms object.
        """

        # Check atype parameter values
        if atype is not None:
            atype = np.asarray(atype, dtype=np.uint64)
            
            # Handle single atype
            if atype.ndim == 0:
                natoms_atype = 1
            
            # Handle list of atype
            elif atype.ndim == 1:
                natoms_atype = atype.shape[0]
            else:
                raise ValueError('Only one atype per atom allowed')
        else:
            atype = np.array([1], dtype=np.uint64)
            natoms_atype = 1
        
        # Check pos parameter values
        if pos is not None:
            pos = np.asarray(pos, dtype=np.float64)
            
            # Handle single pos
            if pos.ndim == 1:
                if pos.shape[0] == 3:
                    natoms_pos = 1
                else:
                    raise ValueError('pos per atom must be 3-dimensional')
            
            # Handle list of pos
            elif pos.ndim == 2:
                if pos.shape[1] == 3:
                    natoms_pos = pos.shape[0]
                else:
                    raise ValueError('pos per atom must be 3-dimensional')
            else:
                raise ValueError('too many dimensions for pos')
        else:
            pos = np.zeros((1, 3), dtype=np.float64)
            natoms_pos = 1
        
        # Check natoms parameter values
        if natoms is not None:
            natoms = int(natoms)
            if ((natoms_atype == 1 or natoms_atype == natoms)
                and (natoms_pos == 1 or natoms_pos == natoms)):
                pass
            else:
                raise ValueError('natoms and length of atype/pos not compatible')
        else:
            if natoms_atype == natoms_pos:
                natoms = natoms_atype
            elif natoms_atype == 1:
                natoms = natoms_pos
            elif natoms_pos == 1:
                natoms = natoms_atype
            else:
                raise ValueError('lengths of atype and pos not compatible')
            
        # Initialize underlying private class attributes
        super(Atoms, self).__setattr__('_Atoms__natoms', natoms)
        super(Atoms, self).__setattr__('_Atoms__prop', Atoms.PropertyDict(self))
        super(Atoms, self).__setattr__('_Atoms__dir', deepcopy(dir(self)))
        
        # Set mandatory properties
        self.prop['atype'] = atype
        self.prop['pos'] = pos
        
        # Set extra properties
        for key, value in iteritems(kwargs):
            self.prop[key] = value
    
    def __str__(self):
        """string output of Atoms data"""
        atype = self.atype
        pos = self.pos
        lines = ['     id |   atype |  pos[0] |  pos[1] |  pos[2]']
        for i in range(self.natoms):
            lines.append('%7i | %7i | %7.3f | %7.3f | %7.3f' % (i, atype[i], pos[i,0], pos[i,1], pos[i,2]))
        
        return '\n'.join(lines)
    
    def __setattr__(self, name, value):
        """Control setting of user-defined attributes"""
        if not hasattr(self, name) or name in self.prop:
            self.prop[name] = value
        else:
            super(Atoms, self).__setattr__(name, value)
    
    def __intslice(self, intnum):
        if intnum == -1:
            return slice(intnum, None)
        else:
            return slice(intnum, intnum+1)
    
    def __deepcopy__(self, memo):
        """Properly handle deepcopy"""
        d = OrderedDict()
        atype = deepcopy(self.prop['atype'])
        pos = deepcopy(self.prop['pos'])
        for key in self.prop:
            if key not in ['atype', 'pos']:
                d[key] = deepcopy(self.prop[key])
        return Atoms(atype=atype, pos=pos, **d)
    
    def __getitem__(self, index):
        """Index getting of Atoms."""
        view = OrderedDict()
        if isinstance(index, inttype):
            index = self.__intslice(index)
            #for key in self.prop.keys():
            #    view[key] = self.prop[key][np.newaxis, index]
        #else:
        for key in self.prop.keys():
            view[key] = self.prop[key][index]
        return Atoms(**view)
        
    def __setitem__(self, index, value):
        """Index setting of Atoms."""
        try:
            assert isinstance(value, Atoms)
            assert sorted(value.prop.keys()) == sorted(self.prop.keys())
        except:
            raise ValueError('Can only set Atoms with matching properties')
        
        if isinstance(index, inttype):
            index = self.__intslice(index)
            
        for key in self.prop.keys():
            self.prop[key][index] = value.prop[key]
    
    def __len__(self):
        """len of atoms = natoms."""
        return self.natoms
    
    @property
    def natoms(self):
        """int : The number of atoms in the Atoms class."""
        return self.__natoms
    
    @property
    def atypes(self):
        """list of int : List of unique atom types."""
        return np.unique(self.atype)
    
    @property
    def natypes(self):
        """int : The number of atom types in the Atoms class."""
        return len(self.atypes)
        
    @property
    def prop(self):
        """PropertyDict : All assigned per-atom properties."""
        return self.__prop
    
    def to_dataframe(self):
        """
        Returns a pandas.DataFrame version of the atoms.  NOTE: The 
        underlying data is copied meaning that calling this at least
        doubles the memory cost.
        """
        
        def indexstr(shape):
            """
            Iterates through all unique indicies of an array with a given shape.
            
            Parameters
            ----------
            shape : tuple of int
                The array shape to iterate through.
                
            Yields
            ------
            index : tuple of int
                A unique index set of the array.
            istr : str
                A string representation of index with numbers in [].
            """
            
            if tuple(shape) == ():
                # Yield for empty shape
                yield (), ''
            else:
                # Loop over all values of the first index of shape
                for i in range(shape[0]):
                    index1 = (i,)
                    istr1 = '['+str(i)+']'
                    
                    # Recursively go through other indicies of shape
                    for index2, istr2 in indexstr(shape[1:]):
                        
                        # Combine index components
                        index = index1 + index2
                        istr = istr1 + istr2
                        
                        # Yield composite
                        yield index, istr
        
        # Initialize new dictionary of values
        values = OrderedDict()
        for key in self.prop.keys():
            
            # Flatten multidimensional arrays
            for index, istr in indexstr(self.prop[key].shape[1:]):
                newkey = key + istr
                
                # Copy values over
                if index == ():
                    values[newkey] = self.prop[key]
                else:
                    values[newkey] = self.prop[key][(Ellipsis, ) + index]
        
        # Return DataFrame
        return pd.DataFrame(values)