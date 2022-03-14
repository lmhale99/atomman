# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from . import Shape, Plane

class Cylinder(Shape):
    """
    Class representing a cylinder in space. 
    """
    def __init__(self,
                 center1: npt.ArrayLike,
                 center2: npt.ArrayLike,
                 radius: float,
                 endcaps: bool = True):
        """
        Defines a cylinder
        
        Parameters
        ----------
        center1 : array-like object
            A point on the cylinder's axis.  If endcaps is True, the point is taken as the center of
            the cylinder's endcap planes.
        center2 : array-like object
            The center point of one of the flat ends.
        radius : float
            The radius
        endcaps : bool, optional
            Indicates if the cylinder is capped at the two
        """
        self.center1 = center1
        self.center2 = center2
        self.radius = radius
        self.endcaps = endcaps
    
    @property 
    def center1(self) -> np.ndarray:
        """numpy.NDArray : 3D vector position of one cylinder end"""
        return self.__center1
    
    @center1.setter
    def center1(self, value: npt.ArrayLike):
        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "center1 must be a 3D vector"
        
        # Save
        self.__center1 = value
    
    @property 
    def center2(self) -> np.ndarray:
        """numpy.NDArray : 3D vector position of one cylinder end"""
        return self.__center2
    
    @center2.setter
    def center2(self, value: npt.ArrayLike):
        
        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "center2 must be a 3D vector"
        
        # Save
        self.__center2 = value
    
    @property
    def radius(self) -> float:
        """float : the cylinder's radius"""
        return self.__radius

    @radius.setter
    def radius(self, value: float):
        # Check that value is proper
        value = float(value)
        assert value > 0, "radius must be positive"

        # Save
        self.__radius = value
    
    @property
    def axis(self) -> np.ndarray:
        """numpy.NDArray : 3D unit vector parallel to cylinder axis."""
        d = self.center2 - self.center1
        return d / np.linalg.norm(d)
    
    @property
    def endcaps(self) -> float:
        """bool : indicates if the endcaps are included in inside/outside determination"""
        return self.__endcaps

    @endcaps.setter
    def endcaps(self, value: bool):
        if not isinstance(value, bool):
            raise ValueError('endcaps must be bool')
        self.__endcaps = value

    def inside(self,
               pos: npt.ArrayLike,
               inclusive: bool = True) -> np.ndarray:
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
        
        axis = self.axis
        
        # Find points <= radius from axis
        distfromaxis = np.linalg.norm(np.cross(pos - self.center1, axis), axis=-1)
        if inclusive:
            insidecircle = distfromaxis <= self.radius
        else:
            insidecircle = distfromaxis < self.radius

        if self.endcaps:
            # Define end planes
            plane1 = Plane(-axis, self.center1)
            plane2 = Plane(axis, self.center2)
        
            # Combine tests
            return insidecircle & plane1.below(pos, inclusive=inclusive) & plane2.below(pos, inclusive=inclusive)
        
        else:
            return insidecircle