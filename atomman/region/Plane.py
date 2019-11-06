# coding: utf-8

# http://www.numpy.org/
import numpy as np

class Plane():
    """
    Class representing a plane in space. Useful for making slices.
    """
    def __init__(self, normal, point):
        """
        Defines a plane in space.
        
        Parameters
        ----------
        normal : array-like object
            3D normal vector of the plane.
        point : array-like object
            3D vector coordinate of any point in the plane.
        """
        
        # Set values
        self.normal = normal 
        self.point = point
    
    @property
    def normal(self):
        """numpy.NDArray : 3D normal unit vector of the plane."""
        return self.__normal
    
    @normal.setter
    def normal(self, value):
        
        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "point must be a 3D vector"
        
        # Save normalvector as a unit vector
        self.__normal = value / np.linalg.norm(value)
    
    @property
    def point(self):
        """numpy.NDArray : a 3D vector position on the plane."""
        return self.__point
    
    @point.setter
    def point(self, value):
        
        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "point must be a 3D vector"
        
        # Save
        self.__point = value
    
    def below(self, pos, inclusive=True):
        """
        Indicates if position(s) are below the plane.  Note that identifying
        points as above or below is dependent on the sign of the plane normal.
        
        Parameters
        ----------
        pos : array-like object
            Nx3 array of coordinates. 
        inclusive : bool, optional
            Indicates if points in the plane are to be included.
            Default value is True.
            
        Returns
        -------
        numpy.NDArray
            N array of bool values: True if below the plane
        """
        
        # Ensure pos is a numpy array 
        pos = np.asarray(pos)
        
        # Determine above/below based on dot product with plane normal vector
        normpoint = np.dot(self.normal, self.point)
        normpos = np.inner(self.normal, pos)
        
        if inclusive:
            return normpos <= normpoint
        else:
            return normpos < normpoint
    
    def above(self, pos, inclusive=False):
        """
        Indicates if position(s) are above the plane.  Note that identifying
        points as above or below is dependent on the sign of the plane normal.
        
        Parameters
        ----------
        pos : array-like object
            Nx3 array of coordinates. 
        inclusive : bool, optional
            Indicates if points in the plane are to be included.
            Default value is False.
        
        Returns
        -------
        numpy.NDArray
            N array of bool values: True if above the plane
        """
        return ~self.below(pos, inclusive=not inclusive)