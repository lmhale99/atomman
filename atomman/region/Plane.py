# coding: utf-8
from __future__ import annotations

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from ..tools import vect_angle

class Plane():
    """
    Class representing a plane in space. Useful for making slices.
    """
    def __init__(self,
                 normal: npt.ArrayLike,
                 point: npt.ArrayLike):
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
    def normal(self) -> np.ndarray:
        """numpy.NDArray : 3D normal unit vector of the plane."""
        return self.__normal

    @normal.setter
    def normal(self, value: npt.ArrayLike):

        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "point must be a 3D vector"

        # Save normalvector as a unit vector
        self.__normal = value / np.linalg.norm(value)

    @property
    def point(self) -> np.ndarray:
        """numpy.NDArray : a 3D vector position on the plane."""
        return self.__point

    @point.setter
    def point(self, value: npt.ArrayLike):

        # Check that value is proper
        value = np.asarray(value)
        assert value.shape == (3,), "point must be a 3D vector"

        # Save
        self.__point = value

    def below(self,
              pos: npt.ArrayLike,
              inclusive: bool = True) -> np.ndarray:
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

    def above(self,
              pos: npt.ArrayLike,
              inclusive: bool = False) -> np.ndarray:
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

    def operate(self,
                rotation: npt.ArrayLike,
                translation: npt.ArrayLike) -> Plane:
        """
        Return a new plane transformed by given symmetry operation

        Parameters
        ----------
        rotation: array-like, (3, 3)
            rotation matrix in cartedian coordinates
        translation: array-like, (3, )
            translation matrix in cartedian coordinates
        
        Returns
        -------
        Plane
            A new plane transformed by the specified operations.
        """
        rotation = np.array(rotation)
        translation = np.array(translation)

        new_normal = np.dot(np.linalg.inv(rotation.T), self.normal)
        new_point = np.dot(rotation, self.point) + translation
        return Plane(new_normal, new_point)

    def __eq__(self, other: Plane) -> bool:
        """
        Compare with a given plane within default tolerance parameters
        """
        if self is other:
            return True

        return self.isclose(other)

    def isclose(self,
                other: Plane,
                atol: float = 1e-8) -> bool:
        """
        Check the plane and a given one represent the same.
        Note that if the normal vectors of the two planes are antiparallel, they are considered to be different.

        Parameters
        ----------
        other: Plane
            a plane to be compared with
        rtol: float
            the relative tolerance parameter used in numpy.isclose
        atol: float
            the absolute tolerance parameter used in numpy.isclose

        Returns
        -------
        bool:
            Return true if the plane and a given one represent the same within tolerance.
        """
        # check if normals are parallel
        angle = vect_angle(self.normal, other.normal, unit='radian')
        if not np.isclose(angle, 0.0, rtol=0.0, atol=atol):
            return False

        # check if the point is contained in the other plane
        projected = np.dot(other.normal, self.point - other.point)
        return np.isclose(projected, 0.0, rtol=0.0, atol=atol)
