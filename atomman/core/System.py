#Standard library imports
from collections import OrderedDict
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
from Atoms import Atoms
from Box import Box
from atomman.tools import nlist, dvect
     
class System:
    #Define atomic system
    
    def __init__(self, atoms=Atoms(), box=Box(), pbc=(True, True, True), scale=False, prop=OrderedDict()):
        #Initialize System instance
                
        assert isinstance(box, Box), 'Invalid box entry'
        self.__box = box
        
        assert isinstance(atoms, Atoms), 'Invalid atoms entry'
        self.__atoms = atoms
        self.atoms_data = self.__atoms.data
        self.atoms_view = self.__atoms.view
        self.atoms_dtype = self.__atoms.dtype
        
        if scale is True:
            self.atoms_view['pos'][:] = self.unscale(atoms.view['pos'])
        
        pbc = np.asarray(pbc, dtype=bool)
        assert pbc.shape == (3L,), 'invalid pbc entry' 
        self.__pbc = pbc
        
        assert isinstance(prop, OrderedDict), 'invalid prop dictionary'
        self.prop = prop

    def atoms(self, arg1=None, arg2=None, arg3=None, dtype=None, unit=None, scale=False): 
    #Returns a copy of atoms if no arguments given, else calls atoms_prop()
    
        if allNone((arg1, arg2, arg3, dtype, unit)):
            if scale:
                newatoms = deepcopy(self.__atoms)
                newatoms.view['pos'][:] = self.scale(newatoms.view['pos'])
                return newatoms
            else:
                return deepcopy(self.__atoms)
        else:
            return self.atoms_prop(arg1=arg1, arg2=arg2, arg3=arg3, dtype=dtype, unit=unit, scale=scale)
        
    def atoms_prop(self, arg1=None, arg2=None, arg3=None, dtype=None, unit=None, scale=False):
    #Calls atoms.prop() method, but with scale option added
  
        #scale=False is identical to atoms.prop()
        if scale is False:
            return self.__atoms.prop(arg1=arg1, arg2=arg2, arg3=arg3, dtype=dtype, unit=unit)
        
        #Determine if a value is given when the arg1 is a prop name
        elif isinstance(arg1, (str, unicode)):
            value = False
            if isinstance(arg2, (int, long)) and arg2 >= 0 and arg2 < self.natoms():
                if arg3 is not None:
                    arg3 = self.unscale(arg3)
                    value = True
            elif arg2 is not None:
                arg2 = self.unscale(arg2)
                value = True
        
        #Determine if a value is given when the arg1 is an atom index
        elif isinstance(arg1, (int, long)) and arg1 >= 0 and arg1 < self.natoms():
            value = False
            if isinstance(arg2, (str, unicode)):
                if arg3 is not None:
                    arg3 = self.unscale(arg3)
                    value = True
            elif arg2 is not None:
                arg2 = self.unscale(arg2)
                value = True
        
        assert unit is None,        'unit and scale cannot be used together'
        
        #Call atoms.prop() using unscaled value
        if value:
            self.__atoms.prop(arg1=arg1, arg2=arg2, arg3=arg3, dtype=dtype)
            self.atoms_data = self.__atoms.data
            self.atoms_view = self.__atoms.view
            self.atoms_dtype = self.__atoms.dtype
        
        #Call atoms.prop() and scale the results
        else:
            return self.scale(self.__atoms.prop(arg1=arg1, arg2=arg2, arg3=arg3, dtype=dtype))              
    
    def natoms(self):
        #Return the number of atoms in the system
        return self.__atoms.natoms()
    
    def natypes(self):
        #Return the max atype value in all of the atoms
        return self.__atoms.natypes()
    
    def box(self, term=None, unit=None):
        #Return a copy of box or call box.get()
        if term is None and unit is None:
            return deepcopy(self.__box)
        else:
            return self.__box.get(term, unit)
    
    def box_get(self, term, unit=None):
        #Access box.get()
        return self.__box.get(term, unit)
    
    def box_set(self, origin=None, vects=None,
                 avect=None, bvect=None, cvect=None,
                 a=None, b=None, c=None, alpha=None, beta=None, gamma=None,
                 lx=None, ly=None, lz=None, xy=None, xz=None, yz=None,
                 xlo=None, xhi=None, ylo=None, yhi=None, zlo=None, zhi=None, unit=None, scale=False): 
        #Set box with scale option
        
        if scale is True:
            spos = self.scale(self.atoms_view['pos'])
        
        self.__box.set(origin=origin, vects=vects,
                       avect=avect, bvect=bvect, cvect=cvect,
                       a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
                       lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz,
                       xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, unit=unit)
        
        if scale is True:
            self.atoms_view['pos'][:] = self.unscale(spos)
    
    def box_normalize(self):
        #Allows for box vectors to be normalized and all atoms adjusted properly
        spos = self.scale(self.atoms_view['pos'])
        self.__box.normalize()
        self.atoms_view['pos'][:] = self.unscale(spos)
    
    def pbc(self, value=None):
        #Get or set pbc conditions
        if value is None:
            return deepcopy(self.__pbc)
        elif isinstance(value, (int, long)):
            return self.__pbc[value]
        else:
            self.__pbc[:] = value
    
    def scale(self, value):
        #Converts from Cartesian to reduced box coordinates

        value = np.asarray(value)
        vects = self.box('vects')
        inverse = np.linalg.inv(vects)
        origin = self.box('origin')
        
        return (value - origin).dot(inverse)
          
    def unscale(self, value):
        #Converts from reduced box to Cartesian coordinates

        value = np.asarray(value)
        vects = self.box('vects')
        origin = self.box('origin')
        
        return value.dot(vects) + origin
    
    def dvect(self, pos_0, pos_1):
        #compute shortest distance between two (sets of) positions taking self system's pbc information into account
        
        pos_0 = np.asarray(pos_0)
        #if pos_0 is integer(s), use as atom indexes of self system 
        if pos_0.dtype == int:
            pos_0 = self.atoms_view['pos'][pos_0]
        
        pos_1 = np.asarray(pos_1)
        #if pos_1 is integer(s), use as atom indexes of self system
        if pos_1.dtype == int:
            pos_1 = self.atoms_view['pos'][pos_1]
        
        #call atomman.tools.dvect using self system's box and pbc
        return dvect(pos_0, pos_1, self.box(), self.pbc())
        
    def wrap(self):
        #Wrap atoms around periodic boundaries and extend non-periodic boundaries if needed
        
        mins = np.array([0.0, 0.0, 0.0])
        maxs = np.array([1.0, 1.0, 1.0])
        spos = self.scale(self.atoms_view['pos'])
        
        for i in xrange(3):
            if self.pbc(i):
                spos[:, i] -= np.floor(spos[:, i])
            else:
                min = spos[:, i].min()
                max = spos[:, i].max()
                if min < 0.0: mins[i] = min - 0.001
                if max > 1.0: maxs[i] = max + 0.001             
                        
        self.atoms_view['pos'][:] = self.unscale(spos)        
        
        origin = self.box('origin') + mins.dot(self.box('vects')) 
        avect = self.box('avect') * (maxs[0] - mins[0])
        bvect = self.box('bvect') * (maxs[1] - mins[1])
        cvect = self.box('cvect') * (maxs[2] - mins[2])
        self.__box = Box(avect=avect, bvect=bvect, cvect=cvect, origin=origin)          
    
    def nlist(self, cutoff, cmult=1):
        #Build neighbor list for all atoms using cutoff.
        self.prop['nlist'] = nlist(self, cutoff, cmult)
                               
def allNone(check):
    for test in check:
        if test is not None:
            return False    
    return True        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
