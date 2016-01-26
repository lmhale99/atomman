#Standard library imports
from collections import OrderedDict
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
import atomman.tools.unitconvert as uc

class Atoms:
    #Class for handling atomic properties
    
    def __init__(self, natoms=1, prop={}, prop_unit=None, prop_dtype=None, nvals=30, data=None, view=None):
        #initialize Atoms instance    
        if data is None and view is None:
            if natoms < 1:
                natoms = 1
            if nvals < 4:
                nvals = 4
            self.data = np.zeros((natoms, nvals), dtype=float)
            
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
                    
                if prop_unit is None:
                    unit = None
                else:
                    unit = prop_unit[key]
                self.prop(key, value, dtype=dtype, unit=unit)
                
        elif data is not None and view is not None and prop_dtype is not None:
            self.data = np.asarray(data, dtype=float)
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
        #string output of Atoms data
        if self.data.ndim == 1:
            atype = np.array([self.view['atype']])
            pos = np.array([self.view['pos']])
        else:
            atype = self.view['atype']
            pos = self.view['pos']
        lines = ['     id |   atype |  pos[0] |  pos[1] |  pos[2]']
        for i in xrange(self.natoms()):
            lines.append('%7i | %7i | %7.3f | %7.3f | %7.3f' % (i, atype[i], pos[i,0], pos[i,1], pos[i,2]))
        
        return '\n'.join(lines)
    
    def __getitem__(self, index):
        view = OrderedDict()
        for k, v in self.view.iteritems():
            view[k] = v[index]
        return Atoms(data=self.data[index], view=view, prop_dtype=self.dtype)
        
    def __setitem__(self, index, value):
        assert isinstance(value, Atoms)
        newdata = self[index].data
        assert value.data.shape == newdata.shape
        assert value.prop() == self.prop()
        assert value.dtype == self.dtype
        newdata[:] = value.data.copy()        
        
    def natoms(self):
        if self.data.ndim == 1:
            return long(1)
        else:
            return long(len(self.data))
        
    def natypes(self):
        return long(self.view['atype'].max())

    def prop(self, arg1=None, arg2=None, arg3=None, dtype=None, unit=None):
    #Accesses the per-atom property lists    
        term = None
        a_id = None
        value = None

        if arg1 is not None:
            if isinstance(arg1, (str, unicode)):
                term = arg1
            elif isinttype(arg1) and arg1 >= 0 and arg1 < self.natoms():
                a_id = arg1
            else:
                raise TypeError('Invalid argument')
        
        if arg2 is not None:
            if isinstance(arg2, (str, unicode)) and term is None:
                term = arg2
            elif isinttype(arg2) and arg2 >= 0 and arg2 < self.natoms() and a_id is None:
                a_id = arg2
            else:
                value = arg2
        if arg3 is not None:
            if value is None:
                value = arg3
            else:
                raise TypeError('Invalid argument')
       
        #Return list of set properties if no term and atom id are given
        if term is None and a_id is None:
            return [k for k in self.view]
        
        #Get or set an atom
        elif term is None:
            if value is None:
                return deepcopy(self[a_id])
            elif isinstance(value, Atoms) and value.natoms() == 1:
                self[a_id] = value
            else:
                raise TypeError('Invalid argument')
                
        #Get or set a property for all atoms
        elif a_id is None:
            if value is None:
                try:
                    if dtype is None:
                        dtype = self.dtype[term]
                    value = np.array(self.view[term], dtype=dtype) 
                except:
                    return None
                
                if value.ndim == 1:
                    if len(value) == 1:
                        value = value[0]
                elif len(value[0]) == 1:
                    value = value.T[0]
                return uc.get_in_units(value, unit)

            else:
                value = np.asarray(value)
                if value.ndim == 1:
                    value = value[:,np.newaxis]
          
                if term not in self.view:
                    assert len(value) == self.natoms(), 'Array size does not match number of atoms'
                    if dtype is None:
                        dtype = value.dtype
                    
                    self.__add_prop(term, value[0].size, value[0].shape, dtype)
                
                #set value if it has the correct dtype 
                assert np.allclose(value, np.asarray(value, dtype=self.dtype[term])), "value not compatible with term's dtype"
                assert dtype is None or dtype==self.dtype[term],                       "stored dtype already assigned to %s" % self.dtype[term]
                
                self.view[term][:] = uc.set_in_units(value, unit)
            
        #Get or set a specific atom's property   
        else:
            if value is None:
                try:
                    if dtype is None:
                        dtype = self.dtype[term]
                    value = np.array(self.view[term][a_id], dtype=dtype) 
                    if len(value) == 1:
                        value = value[0]
                    return uc.get_in_units(value, unit)
                except:
                    return None
            else:
                value = np.asarray(value)
          
                if term not in self.view:
                    if dtype is None:
                        dtype = value.dtype
                    
                    self.__add_prop(term, value.size, value.shape, dtype)
                
                #set value if it has the correct shape and dtype 
                assert np.allclose(value, np.asarray(value, dtype=self.dtype[term])), "value not compatible with term's dtype"
                assert dtype is None or dtype==self.dtype[term],                       "stored dtype already assigned to %s" % self.dtype[term]

                self.view[term][a_id] = uc.set_in_units(value, unit)
                
    def __add_prop(self, term, size, shape, dtype):
        try:
            end = self.__start + size
        except:
            raise TypeError('New properties can only be added to ')
                
        #check if self.data is too small
        if end > len(self.data[0]):
            
            #create larger array
            data = np.empty((self.natoms(), end +  10))
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
        self.view[term] = self.data[:,  self.__start : end]
        self.view[term].shape = (self.natoms(), ) + shape
        self.dtype[term] = dtype
        self.__start = end
                
                
def isinttype(value):
    if isinstance(value, (int, long)): 
        return True
    elif isinstance(value, np.ndarray) and len(value) == 0:
        if value.dtype == int or value.dtype == long: 
            return True
        else:
            try:
                if issubclass(value.dtype.type, np.int):
                    return True
                else:
                    return False
            except:
                return False
    else:
        return False                 
                
                
                
                
                
                
                
                
                
                
                
                
            
        
        