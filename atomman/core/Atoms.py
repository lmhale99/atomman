#Standard library imports
from collections import OrderedDict
from copy import deepcopy

#External library imports
import numpy as np

from atomman.tools import is_dtype_int, is_dtype_bool

class Atoms(object):
    """
    Class for handling atomic properties.
    
    Attributes:
    natoms -- number of atoms.
    natypes -- number of assigned atom types.
    data -- the underlying numpy array used.
    view -- dictionary of structured views for each assigned property
    dtype -- dictionary of data types for each assigned property
    
    Methods:
    prop -- accesses the atom properties for value assignment or retrieval.
    """
    
    def __init__(self, natoms=1, prop={}, prop_dtype=None, nvals=30, data=None, view=None):
        """
        Initilizes an Atoms instance.
        
        Keyword Arguments:
        natoms -- Number of atoms.  This is fixed once the system is initilized.  Default value is 1.
        prop -- Dictionary containing property terms.  Allows property values to be defined during initilization.
        prop_dtype -- Dictionary that allows for the data type of properties to be explicitly assigned.
        nvals -- Inititial number of values per atom. Will expand if necessary. Default value is 30.
        data -- numpy array containing all the property data for all the atoms. This is for more explicit initilization control.
        view -- Dictionary containing views of data relating to each property. This is for more explicit initilization control. 
        """    
        if data is None and view is None:
            if natoms < 1:
                natoms = 1
            if nvals < 4:
                nvals = 4
            self.data = np.zeros((natoms, nvals), dtype='float64')
            
            self.view  = OrderedDict()
            self.dtype = OrderedDict()
        
            self.view['atype'] =  self.data[:, 0:1]
            self.view['atype'][:] = np.zeros((natoms, 1))
            self.dtype['atype'] = np.dtype('int32')
        
            self.view['pos'] =  self.data[:, 1:4]
            self.view['pos'][:] = np.zeros((natoms, 3))
            self.dtype['pos'] = np.dtype('float64')
        
            self.__start = 4
        
            for key, value in prop.iteritems():
                if prop_dtype is None:
                    dtype = None
                else:
                    dtype = prop_dtype[key]
                    
                self.prop(key=key, value=value, dtype=dtype)
                
        elif data is not None and view is not None and prop_dtype is not None:
            self.data = np.asarray(data, dtype='float64')
            assert isinstance(view, OrderedDict),       'view must be an OrderedDict'
            assert isinstance(prop_dtype, OrderedDict), 'prop_dtype must be an OrderedDict'
            assert 'atype' in view,                     'view does not have atype'
            assert 'pos' in view,                       'view does not have pos'
            for k in view:
                assert k in prop_dtype,                 'view and prop_dtype keys do not match'
            for k in prop_dtype:
                assert k in view,                       'view and prop_dtype keys do not match'
            
            self.view = view
            self.dtype = prop_dtype  
            self.__start = 0
            for k, v in view.iteritems():
                self.__start += v[0].size
        
    def __str__(self):
        """string output of Atoms data"""
        if self.data.ndim == 1:
            atype = np.array([self.view['atype']])
            pos = np.array([self.view['pos']])
        else:
            atype = self.view['atype']
            pos = self.view['pos']
        lines = ['     id |   atype |  pos[0] |  pos[1] |  pos[2]']
        for i in xrange(self.natoms):
            lines.append('%7i | %7i | %7.3f | %7.3f | %7.3f' % (i, atype[i], pos[i,0], pos[i,1], pos[i,2]))
        
        return '\n'.join(lines)
    
    def __getitem__(self, index):
        """Index getting of atoms."""
        view = OrderedDict()
        for k, v in self.view.iteritems():
            view[k] = v[index]
        return Atoms(data=self.data[index], view=view, prop_dtype=self.dtype)
        
    def __setitem__(self, index, value):
        """Index setting of atoms."""
        assert isinstance(value, Atoms)
        assert value.data.shape == self[index].data.shape
        assert value.prop() == self.prop()
        assert value.dtype == self.dtype
        self.data[index] = value.data.copy()       
    
    @property
    def natoms(self):
        """The number of atoms."""
        if self.data.ndim == 1:
            return long(1)
        else:
            return long(len(self.data))
    
    @property
    def natypes(self):
        """The maximum atype value."""
        return long(self.view['atype'].max())

    def prop(self, **kwargs):
        """
        Accesses the per-atom properties for getting and setting.
        
        Keyword Arguments:
        a_id -- atom index.
        key -- atom property name.
        value -- value(s) to assign to properties associated with a_id and/or term.
        dtype -- data type to explicitly set or retrieve value as. 

        If no arguments given, returns a list of the assigned property keys. Otherwise, a_id and/or key must be specified. The key specifies which property, and the a_id which atom(s) to access. With no value argument, prop() returns which value(s) are associated with the given a_id and/or key. With a value argument, the value is saved according to the given a_id and/or key.
        """                  
        key = kwargs.pop('key', None)
        a_id = kwargs.pop('a_id', None)
        value = kwargs.pop('value', None)
        dtype = kwargs.pop('dtype', None)
        assert len(kwargs) == 0, 'Invalid keyword arguments'
        
        #Return list of assigned properties when no arguments given.
        if key is None and a_id is None:
            assert dtype is None, 'dtype cannot be specified without key'
            assert value is None, 'value cannot be specified without key and/or a_id'
            return [k for k in self.view]
       
        #Get or set an atom
        if key is None and a_id is not None:
            assert dtype is None, 'dtype cannot be specified without key'
            if value is None:
                return deepcopy(self[a_id])
            elif isinstance(value, Atoms) and value.natoms == 1:
                self[a_id] = value
            else:
                raise TypeError('value for a_id must be an Atoms instance of length 1')
                
        #Get or set a property for all atoms
        elif a_id is None and key is not None:
            if value is None:
                try:
                    if dtype is None:
                        dtype = self.dtype[key]
                    value = np.array(self.view[key], dtype=dtype) 
                except:
                    return None
                
                if value.ndim == 1:
                    if len(value) == 1:
                        value = value[0]
                elif len(value[0]) == 1:
                    value = value.T[0]
                return value

            else:
                value = np.asarray(value)

                if value.ndim == 0:
                    value = value.reshape(1)
                if value.ndim == 1:
                    value = value[:,np.newaxis]

                if key not in self.view:
                    if dtype is None:
                        dtype = value.dtype
                    #assert len(value) == self.natoms, 'Array size does not match number of atoms'                    
                    self.__add_prop(key, value[0].size, value[0].shape, dtype)
                if dtype is not None:
                    assert dtype == self.dtype[key], 'dtype already assigned as %s' % self.dtype[key]
                
                #set value if it has the correct dtype 
                if is_dtype_int(self.dtype[key]): 
                    assert is_dtype_int(value.dtype) or is_dtype_bool(value.dtype), 'value is incompatible with int data type'
                elif is_dtype_bool(self.dtype[key]): 
                    assert is_dtype_bool(value.dtype), 'value is incompatible with bool data type'
                self.view[key][:] = value
            
        #Get or set a specific atom's property   
        else:
            if value is None:
                try:
                    if dtype is None:
                        dtype = self.dtype[key]
                    value = np.array(self.view[key], dtype=dtype) 
                except:
                    return None
                value = value[a_id]
                if len(value) == 1:
                    value = value[0]
                return value
            else:
                value = np.asarray(value)
                if value.ndim == 0:
                    value = value.reshape(1)
                
                if key not in self.view:
                    if dtype is None:
                        dtype = value.dtype
                    self.__add_prop(key, value.size, value.shape, dtype)
                if dtype is not None:
                    assert dtype == self.dtype[key], 'dtype already assigned as %s' % self.dtype[key]
                
                #set value if it has the correct shape and dtype 
                if is_dtype_int(self.dtype[key]): 
                    assert is_dtype_int(value.dtype) or is_dtype_bool(value.dtype), 'value is incompatible with int data type'
                elif is_dtype_bool(self.dtype[key]): 
                    assert is_dtype_bool(value.dtype), 'value is incompatible with bool data type'

                self.view[key][a_id] = value
                
    def __add_prop(self, key, size, shape, dtype):
        """Internal function that adds a new property and extends underlying array if needed"""
        try:
            end = self.__start + size
        except:
            raise TypeError('New properties can only be added to ')
                
        #check if self.data is too small
        if end > len(self.data[0]):
            
            #create larger array
            data = np.zeros((self.natoms, end +  10))
            data[:, :self.__start] = self.data[:, :self.__start]
            self.data = data
            
            #reassign views to new array
            start = 0
            for k in self.view:
                vshape = self.view[k].shape
                self.view[k] = self.data[:, start : start + self.view[k][0].size]
                self.view[k].shape = vshape
                start = start + self.view[k][0].size
            assert start == self.__start, 'Error reassigning properties'
            
        #create view and save dtype
        self.view[key] = self.data[:,  self.__start : end]
        self.view[key].shape = (self.natoms, ) + shape
        self.dtype[key] = dtype
        self.__start = end