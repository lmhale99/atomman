#Standard library imports
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
from atomman.tools import vect_angle

class Box(object):
    """
    Class for representing a triclinic box. All property attributes can be retrieved individually, but must be set using one of the given set* methods.
    
    Property Attributes:
    vects -- 3x3 array of box vectors. Can be set directly.
    avect, bvect, cvect -- individual box vectors.
    a, b, c, alpha, beta, gamma -- crystal cell parameter representation of a box.
    lx, ly, lz -- LAMMPS-style box lengths.
    xy, xz, yz -- LAMMPS-style box tilt factors.
    xlo, xhi, ylo, yhi, zlo, zhi -- LAMMPS-style box hi/los.
    origin -- box origin from which the vects are added to define the box. Can be set directly.
    
    Methods:
    set() -- calls one of the other set_* methods if keyword arguments form a full parameter set.
    set_vectors() -- set box with avect, bvect, cvect, (and origin).
    set_abc() -- set box with a, b, c, (alpha, beta, gamma, and origin).
    set_lengths() -- set box with lx, ly, lz, (xy, xz, yz, and origin).
    set_hi_los() -- set box with xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).
    normalize() -- normalizes a generic box to be compatible with the LAMMPS representations.
    """

    def __init__(self, **kwargs):
        """
        Initilizes a box.
        
        If no arguments given, box is set to square unit box with origin = [0,0,0].
        Otherwise, must give one of the following sets:
        - vects, (and origin).
        - avect, bvect, cvect, (and origin).
        - a, b, c, (alpha, beta, gamma, and origin).
        - lx, ly, lz, (xy, xz, yz, and origin).
        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).
        """
        self.__vects = np.eye(3, dtype='float64')
        self.__origin = np.zeros(3, dtype='float64')
        self.__isnorm = True
        
        if len(kwargs) > 0:
            self.set(**kwargs)
    
    @property
    def vects(self):
        """3x3 array of box vectors."""
        return deepcopy(self.__vects)
        
    @vects.setter
    def vects(self, value):
        
        self.__vects[:] = value
        #Zero out near zero terms
        self.__vects[np.isclose(self.__vects/self.__vects.max(), 0.0, atol=1e-9)] = 0.0
        
        #Check if box vectors are consistent with LAMMPS normalization
        if (self.__vects[0,1] == 0.0 and self.__vects[0,2] == 0.0 and self.__vects[1,2] == 0.0 and
            self.__vects[0,0] > 0.0 and self.__vects[1,1] > 0.0 and self.__vects[2,2] > 0.0 and
            abs(self.__vects[1,0]) < self.__vects[0,0]/2. and 
            abs(self.__vects[2,0]) < self.__vects[0,0]/2. and 
            abs(self.__vects[2,1]) < self.__vects[1,1]/2.):
            self.__isnorm = True
        else:
            self.__isnorm = False
        
    @property
    def origin(self):
        """box origin vector position"""
        return deepcopy(self.__origin)
        
    @origin.setter
    def origin(self, value):
        self.__origin[:] = value
        
    @property
    def avect(self):
        """box a vector"""
        return self.vects[0]
        
    @property
    def bvect(self):
        """box b vector"""
        return self.vects[1]
    
    @property
    def cvect(self):
        """box c vector"""
        return self.vects[2]
    
    @property
    def a(self):
        """cell length a, i.e. magnitude of box a vector"""
        return (self.__vects[0,0]**2 + self.__vects[0,1]**2 + self.__vects[0,2]**2)**0.5

    @property
    def b(self):
        """cell length b, i.e. magnitude of box b vector"""
        return (self.__vects[1,0]**2 + self.__vects[1,1]**2 + self.__vects[1,2]**2)**0.5       
    
    @property
    def c(self):
        """cell length c, i.e. magnitude of box c vector"""
        return (self.__vects[2,0]**2 + self.__vects[2,1]**2 + self.__vects[2,2]**2)**0.5 
    
    @property    
    def alpha(self):
        """cell angle alpha in degrees"""
        return vect_angle(self.__vects[1], self.__vects[2])
        
    @property    
    def beta(self):
        """cell angle beta in degrees"""
        return vect_angle(self.__vects[0], self.__vects[2])

    @property    
    def gamma(self):
        """cell angle gamma in degrees"""
        return vect_angle(self.__vects[0], self.__vects[1])      

    @property    
    def lx(self):
        """LAMMPS-style box length lx"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[0,0]
    
    @property    
    def ly(self):
        """LAMMPS-style box length ly"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,1]

    @property    
    def lz(self):
        """LAMMPS-style box length lz"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,2]
        
    @property    
    def xy(self):
        """LAMMPS-style box tilt factor xy"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,0]
        
    @property    
    def xz(self):
        """LAMMPS-style box tilt factor xz"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,0]

    @property    
    def yz(self):
        """LAMMPS-style box tilt factor yz"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,1]        
         
    @property    
    def xlo(self):
        """LAMMPS-style box xlo"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0]
        
    @property    
    def ylo(self):
        """LAMMPS-style box ylo"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1]

    @property    
    def zlo(self):
        """LAMMPS-style box zlo"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2]   
        
    @property    
    def xhi(self):
        """LAMMPS-style box xhi"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0] + self.__vects[0,0]
        
    @property    
    def yhi(self):
        """LAMMPS-style box yhi"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1] + self.__vects[1,1]

    @property    
    def zhi(self):
        """LAMMPS-style box zhi"""
        assert self.__isnorm,                       'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2] + self.__vects[2,2]         
    
    def __str__(self):
        """Returns a string representation of the box"""
        return '\n'.join(['avect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[0,0], self.__vects[0,1], self.__vects[0,2]),
                          'bvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[1,0], self.__vects[1,1], self.__vects[1,2]),
                          'cvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[2,0], self.__vects[2,1], self.__vects[2,2]),
                          'origin = [%6.3f, %6.3f, %6.3f]' % (self.__origin[0],  self.__origin[1],  self.__origin[2])])
    
    
    def set(self, **kwargs):            
        """
        Set the box values.
        
        If no arguments given, box is set to square unit box with origin = [0,0,0].
        Otherwise, must give one of the following sets:
        - vects, (and origin).
        - avect, bvect, cvect, (and origin).
        - a, b, c, (alpha, beta, gamma, and origin).
        - lx, ly, lz, (xy, xz, yz, and origin).
        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).
        """
                    
        #Set default values if no kwargs given
        if len(kwargs) == 0:
            self.vects = np.eye(3)
            self.origin = np.zeros(3)
        
        #Set directly if vects given
        elif 'vects' in kwargs:
            vects = kwargs.pop('vects')
            origin = kwargs.pop('origin', [0.0, 0.0, 0.0])
            assert len(kwargs) == 0, 'Invalid arguments'
            self.vects = vects
            self.origin = origin
        
        #Call set_vectors if vect inputs given
        elif 'avect' in kwargs:            
            self.set_vectors(**kwargs)
        
        #Call set_lengths if length inputs given
        elif 'lx' in kwargs:           
            self.set_lengths(**kwargs)
    
        #Call set_hi_los if hi/lo inputs given
        elif 'xlo' in kwargs:            
            self.set_hi_los(**kwargs)
        
        #Call set_abc if vector magnitudes are given
        elif 'a' in kwargs:            
            self.set_abc(**kwargs)
        
        else:
            raise TypeError('Invalid arguments')
        
    def set_vectors(self, avect, bvect, cvect, origin=[0.0, 0.0, 0.0]):
        """Set the box using vectors (avect, bvect, cvect) and origin."""
        
        self.vects = [avect, bvect, cvect]
        self.origin = origin
        
    def set_abc(self, a, b, c, alpha=90.0, beta=90.0, gamma=90.0, origin=[0.0, 0.0, 0.0]):
        """Set the box using crystal-style cell parameters (a, b, c, alpha, beta, gamma) and origin."""
         
        #Convert to lx, ly, lz, xy, xz, yz
        lx = a
        xy = b * np.cos(gamma * np.pi / 180)
        xz = c * np.cos(beta * np.pi / 180)
        ly = (b**2 - xy**2)**0.5
        yz = (b * c * np.cos(alpha * np.pi / 180) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5

        #Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)
        
    def set_lengths(self, lx, ly, lz, xy=0.0, xz=0.0, yz=0.0, origin=[0.0, 0.0, 0.0]):
        """Set the box using LAMMPS-style lengths (lx, ly, lz), tilt factors (xy, xz, yz) and origin."""
        
        assert lx > 0 and ly > 0 and lz > 0,                            'All lengths must be positive'
        assert abs(xy) < lx/2. and abs(xz) < lx/2. and abs(yz) < ly/2., 'Invalid angles/tilts'
        
        #Construct vects array
        self.vects = [[lx, 0.0, 0.0],
                      [xy, ly,  0.0],
                      [xz, yz,  lz]]
        self.origin = origin
    
    def set_hi_los(self, xlo, xhi, ylo, yhi, zlo, zhi, xy=0.0, xz=0.0, yz=0.0):
        """Set the box using LAMMPS-style hi/los (xlo, xhi, ylo, yhi, zlo, zhi) and tilt factors (xy, xz, yz)."""   

        #Convert to hi/los to lengths and origin
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        origin = [xlo, ylo, zlo]
        
        #Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)
    
    def normalize(self):
        """Rotates the box so that avect_y = avect_z = bvect_z = 0.0."""
        
        #Extract vectors and vector magnitudes
        avect = self.avect
        bvect = self.bvect
        cvect = self.cvect
        a = (avect[0]**2 + avect[1]**2 + avect[2]**2)**0.5 
        b = (bvect[0]**2 + bvect[1]**2 + bvect[2]**2)**0.5 
        c = (cvect[0]**2 + cvect[1]**2 + cvect[2]**2)**0.5 
        
        #test right handedness
        test = np.dot(np.cross(avect, bvect), cvect)
        assert test > 0, 'Supplied vectors are not right handed'
        
        #Convert to normalized lx, ly, lz, xy, xz, yz
        lx = a
        xy = np.dot(bvect, avect/lx)
        ly = (b**2 - xy**2)**0.5
        xz = np.dot(cvect, avect/lx)
        yz = (np.dot(bvect, cvect) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5
    
        #Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=self.origin)     