# coding: utf-8

# http://www.numpy.org/
import numpy as np

from . import Shape

class Sphere(Shape):
    """
    Class representing a sphere in space.
    """
    def __init__(self, center, radius):
        """
        Defines a sphere

        Parameters
        ----------
        center : array-like object
            The position of the sphere's center.
        radius : float
            The sphere's radius. 
        """
        self.center = center
        self.radius = radius
    
    @property 
    def center(self):
        """numpy.NDArray : 3D vector position of sphere's center"""
        return self.__center
    
    @center.setter
    def center(self, value):
        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "center must be a 3D vector"
        
        # Save
        self.__center = value
    
    @property
    def radius(self):
        """float : the sphere's radius"""
        return self.__radius

    @radius.setter
    def radius(self, value):
        # Check that value is proper
        value = float(value)
        assert value > 0, "radius must be positive"

        # Save
        self.__radius = value
        
    def inside(self, pos, inclusive=True):
        """
        Indicates if position(s) are inside the shape.
        
        Parameters
        ----------
        pos : array-like object
            Nx3 array of coordinates. 
        inclusive : bool, optional
            Indicates if points on the shape's boundaries are to be included.
            Default value is True.
        
        Returns
        -------
        numpy.NDArray
            N array of bool values: True if inside shape
        """
        
        # Ensure pos is a numpy array 
        pos = np.asarray(pos)
        
        # Find distance of pos from center
        rpos = np.linalg.norm(pos - self.center, axis=-1)
        
        if inclusive:
            return rpos <= self.radius
        else:
            return rpos < self.radius