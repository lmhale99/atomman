#Standard library imports
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
from atomman.tools import vect_angle
import atomman.tools.unitconvert as uc

class Box:
#Class for representing a triclinic box

    def __init__(self, origin=None, vects=None,
                 avect=None, bvect=None, cvect=None,
                 a=None, b=None, c=None, alpha=None, beta=None, gamma=None,
                 lx=None, ly=None, lz=None, xy=None, xz=None, yz=None,
                 xlo=None, xhi=None, ylo=None, yhi=None, zlo=None, zhi=None, unit=None):
        
        self.__vects = np.eye(3)
        self.__origin = np.zeros(3)
        self.__isnorm = True
       
        self.set(origin=origin, vects=vects,
                 avect=avect, bvect=bvect, cvect=cvect,
                 a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
                 lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz,
                 xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, unit=unit)
    
    def __str__(self):
        return '\n'.join(['avect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[0,0], self.__vects[0,1], self.__vects[0,2]),
                          'bvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[1,0], self.__vects[1,1], self.__vects[1,2]),
                          'cvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[2,0], self.__vects[2,1], self.__vects[2,2]),
                          'origin = [%6.3f, %6.3f, %6.3f]' % (self.__origin[0],  self.__origin[1],  self.__origin[2])])
    
    
    def set(self, origin=None, vects=None,
                 avect=None, bvect=None, cvect=None,
                 a=None, b=None, c=None, alpha=None, beta=None, gamma=None,
                 lx=None, ly=None, lz=None, xy=None, xz=None, yz=None,
                 xlo=None, xhi=None, ylo=None, yhi=None, zlo=None, zhi=None, unit=None):            
        
        #Call set_vectors if vect inputs given
        if vects is not None or allnotNone((avect, bvect, cvect)):
            assert allNone((lx, ly, lz, xy, xz, yz)),                       'Invalid parameter set' 
            assert allNone((xlo, xhi, ylo, yhi, zlo, zhi)),                 'Invalid parameter set' 
            assert allNone((a, b, c, alpha, beta, gamma)),                  'Invalid parameter set'
            
            self.set_vectors(origin=origin, vects=vects, avect=avect, bvect=bvect, cvect=cvect, unit=unit)
        
        #Call set_lengths if length inputs given
        elif allnotNone((lx, ly, lz)):
            assert allNone((vects, avect, bvect, cvect)),                   'Invalid parameter set'
            assert allNone((xlo, xhi, ylo, yhi, zlo, zhi)),                 'Invalid parameter set' 
            assert allNone((a, b, c, alpha, beta, gamma)),                  'Invalid parameter set'
           
            self.set_lengths(origin=origin, lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, unit=unit)
    
        #Call set_hi_los if hi/lo inputs given
        elif allnotNone((xlo, xhi, ylo, yhi, zlo, zhi)):
            assert allNone((vects, avect, bvect, cvect)),                   'Invalid parameter set'
            assert allNone((lx, ly, lz)),                                   'Invalid parameter set' 
            assert allNone((a, b, c, alpha, beta, gamma, origin)),          'Invalid parameter set'
            
            self.set_hi_los(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, xy=xy, xz=xz, yz=yz, unit=unit)
        
        #Call set_abc if vector magnitudes are given
        elif allnotNone((a, b, c)):
            assert allNone((vects, avect, bvect, cvect)),                   'Invalid parameter set'
            assert allNone((lx, ly, lz, xy, xz, yz)),                       'Invalid parameter set'         
            assert allNone((xlo, xhi, ylo, yhi, zlo, zhi)),                 'Invalid parameter set' 
            
            self.set_abc(origin=origin, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, unit=unit)
            
        #Set to cubed unit box if no parameters given    
        else:
            assert allNone((vects, avect, bvect, cvect)),                   'Invalid parameter set'
            assert allNone((lx, ly, lz, xy, xz, yz)),                       'Invalid parameter set'         
            assert allNone((xlo, xhi, ylo, yhi, zlo, zhi)),                 'Invalid parameter set' 
            assert allNone((a, b, c, alpha, beta, gamma, origin)),          'Invalid parameter set'
            
            self.__vects[:] = uc.set_in_units(np.eye(3), unit)
            self.__origin[:] = np.zeros(3)
            self.__isnorm = True
    
    def set_origin(self, origin=None, unit=None):
    #Sets the box origin    
        
        #Default origin is (0,0,0)
        if origin is None:
            origin = np.zeros(3)
        
        #Check shape of origin
        origin = np.asarray(origin)
        assert origin.shape == (3L,)
        
        #Convert units and save
        self.__origin[:] = uc.set_in_units(origin, unit)
        
    def set_vectors(self, origin=None, vects=None, avect=None, bvect=None, cvect=None, unit=None):
    #Set the direction vectors of the box. 
        
        #Test if vects or (avect, bvect, cvect) are given and build vects
        if vects is not None:
            assert allNone((avect, bvect, cvect)),                'cannot specify vects with avect, bvect, or cvect'
            vects = np.asarray(vects)
        else:
            assert allnotNone((avect, bvect, cvect)),             'all avect, bvect, cvect needed'
            vects = np.array([avect, bvect, cvect])
        
        #Check shape of vects
        assert vects.shape == (3L, 3L)
        
        #Zero out near zero terms
        vects[np.isclose(vects/vects.max(), 0.0)] = 0.0
        
        #Convert units and save
        self.__vects[:] = uc.set_in_units(vects, unit)
        
        #Check if box vectors are consistent with LAMMPS normalization
        if (vects[0,1] == 0.0 and vects[0,2] == 0.0 and vects[1,2] == 0 and
            vects[0,0] > 0.0 and vects[1,1] > 0.0 and vects[2,2] > 0.0 and
            abs(vects[1,0]) < vects[0,0]/2. and 
            abs(vects[2,0]) < vects[0,0]/2. and 
            abs(vects[2,1]) < vects[1,1]/2.):
            self.__isnorm = True
        else:
            self.__isnorm = False
        
        #Call set_origin
        self.set_origin(origin, unit)
        
    def set_abc(self, origin=None, a=None, b=None, c=None, alpha=None, beta=None, gamma=None, unit=None):
    #Set the direction vectors using vector magnitudes (a, b, c) and angles (alpha, beta, gamma). 
    
        #Unit box constructed if no a,b,c are given
        if allNone((a, b, c)):
            a = b = c = 1.0
        assert allnotNone((a, b, c)),   'All or none a,b,c must be given' 

        #Default angles are 90.0 degrees
        if alpha is None: alpha = 90.0
        if beta is None:  beta =  90.0
        if gamma is None: gamma = 90.0
        
        #Convert to lx, ly, lz, xy, xz, yz
        lx = a
        xy = b * np.cos(gamma * np.pi / 180)
        xz = c * np.cos(beta * np.pi / 180)
        ly = (b**2 - xy**2)**0.5
        yz = (b * c * np.cos(alpha * np.pi / 180) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5

        #Call set_lengths
        self.set_lengths(origin=origin, lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, unit=unit)
        
    def set_lengths(self, origin=None, lx=None, ly=None, lz=None, xy=None, xz=None, yz=None, unit=None):
    #Set the direction vectors using lengths (lx, ly, lz) and tilts (xy, xz, yz).
        
        #Check that all lengths are given and positive
        assert allnotNone((lx, ly, lz)),                                'All lengths lx, ly, and lz must be given' 
        assert lx > 0 and ly > 0 and lz > 0,                            'All lengths must be positive'
        
        #Default tilts are 0.0
        if xy is None: xy = 0.0
        if xz is None: xz = 0.0
        if yz is None: yz = 0.0
        assert abs(xy) < lx/2. and abs(xz) < lx/2. and abs(yz) < ly/2., 'Invalid angles/tilts'
        
        #Construct vects array
        vects = np.array([[lx, 0.0, 0.0],
                          [xy, ly,  0.0],
                          [xz, yz,  lz]], dtype=float)
        vects[np.isclose(vects/vects.max(), 0.0)] = 0.0
        
        #Save in appropriate units
        self.__vects[:] = uc.set_in_units(vects, unit)

        #vects is constructed to be LAMMPS normalized
        self.__isnorm = True
        
        #Call set_origin
        self.set_origin(origin, unit)
    
    def set_hi_los(self, xlo=None, xhi=None, ylo=None, yhi=None, zlo=None, zhi=None, xy=None, xz=None, yz=None, unit=None):
    #Set the direction vectors using hi/los (xlo, xhi, ylo, yhi, zlo, zhi) and tilts (xy, xz, yz).    
        
        #Check that all hi/los are given
        assert allnotNone((xlo, xhi, ylo, yhi, zlo, zhi)),        'All hi/los must be given'
        
        #Convert to lengths
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        
        #Build origin
        origin = np.array([xlo, ylo, zlo], dtype=float)
        
        #Call set_lengths
        self.set_lengths(origin=origin, lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, unit=unit)
    
    def normalize(self):
    #Rotates the box so that avect_y = avect_z = bvect_z = 0.
        
        #test right handedness
        test = np.dot(np.cross(self.__vects[0], self.__vects[1]), self.__vects[2])
        assert test > 0, 'Supplied vectors are not right handed'
        
        #Extract vectors and vector magnitudes
        avect = self.__vects[0]
        bvect = self.__vects[1]
        cvect = self.__vects[2]
        a = (avect[0]**2 + avect[1]**2 + avect[2]**2)**0.5 
        b = (bvect[0]**2 + bvect[1]**2 + bvect[2]**2)**0.5 
        c = (cvect[0]**2 + cvect[1]**2 + cvect[2]**2)**0.5 
        
        #Convert to normalized lx, ly, lz, xy, xz, yz
        lx = a
        xy = np.dot(bvect, avect/lx)
        ly = (b**2 - xy**2)**0.5
        xz = np.dot(cvect, avect/lx)
        yz = (np.dot(bvect, cvect) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5
    
        #Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=self.__origin)
            
    def get(self, term, unit=None):
    #Returns vectors, magnitudes, angles, lengths, hi/los and tilts
        
        #vects, avect, bvect, cvect
        if term == 'vects':
            value = deepcopy(self.__vects)
        elif term == 'avect':
            value = deepcopy(self.__vects[0])
        elif term == 'bvect':
            value = deepcopy(self.__vects[1])
        elif term == 'cvect':
            value = deepcopy(self.__vects[2])
            
        #a, b, c, alpha, beta, gamma    
        elif term == 'a':
            value = (self.__vects[0,0]**2 + self.__vects[0,1]**2 + self.__vects[0,2]**2)**0.5
        elif term == 'b':
            value = (self.__vects[1,0]**2 + self.__vects[1,1]**2 + self.__vects[1,2]**2)**0.5
        elif term == 'c':
            value = (self.__vects[2,0]**2 + self.__vects[2,1]**2 + self.__vects[2,2]**2)**0.5
        elif term == 'alpha':
            value = vect_angle(self.__vects[1], self.__vects[2])
        elif term == 'beta':
            value = vect_angle(self.__vects[0], self.__vects[2])
        elif term == 'gamma':
            value = vect_angle(self.__vects[0], self.__vects[1])
            
        #lx, ly, lz    
        elif term == 'lx':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[0,0]
        elif term == 'ly':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[1,1]
        elif term == 'lz':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[2,2]
            
        #xy, xz, yz    
        elif term == 'xy':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[1,0]
        elif term == 'xz':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[2,0]
        elif term == 'yz':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__vects[2,1]
            
        #xlo, ylo, zlo    
        elif term == 'xlo':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[0]
        elif term == 'ylo':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[1]    
        elif term == 'zlo':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[2]
            
        #xhi, yhi, zhi
        elif term == 'xhi':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[0] + self.__vects[0,0]
        elif term == 'yhi':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[1] + self.__vects[1,1]
        elif term == 'zhi':
            assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
            value = self.__origin[2] + self.__vects[2,2]  
        
        #origin
        elif term == 'origin':
            value = deepcopy(self.__origin)
        
        #unknown
        else:
            return None
            
        return uc.get_in_units(value, unit)

def allNone(check):
    for test in check:
        if test is not None:
            return False    
    return True

def allnotNone(check):
    for test in check:
        if test is None:
            return False
    return True
    