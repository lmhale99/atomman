# coding: utf-8

# Standard Python libraries
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from .. import Box, ElasticConstants
from ..tools import axes_check, vect_angle, miller

class VolterraDislocation(object):
    """
    Generic class for a Volterra solution of a straight dislocation.
    """

    def __init__(self,
                 C: ElasticConstants,
                 burgers: npt.ArrayLike,
                 ξ_uvw: Optional[npt.ArrayLike] = None,
                 slip_hkl: Optional[npt.ArrayLike] = None,
                 transform: Optional[npt.ArrayLike] = None,
                 axes: Optional[npt.ArrayLike] = None,
                 box: Optional[Box] = None,
                 m: Union[str, npt.ArrayLike] = 'x',
                 n: Union[str, npt.ArrayLike] = 'y',
                 cart_axes: bool = False,
                 tol: float = 1e-8):
        """
        Initializes the solution. Calls solve.

        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Burgers vector.
        ξ_uvw : array-like object
            The Miller crystal vector associated with the dislocation's line
            direction.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        slip_hkl : array-like object
            The Miller plane indices associated with the dislocation's slip
            plane.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        transform : array-like object, optional
            A 3x3 set of orthogonal Cartesian vectors that define the
            transformation matrix to use on C and burgers to convert from the
            standard (unit cell) orientation to the dislocation orientation.
            The 3 vectors will automatically be converted into unit vectors.
            Using this is an alternative to using ξ_uvw and slip_hkl.
        axes : array-like object, optional
            Same as transform.  Retained for backwards compatibility.
        box : atomman.Box, optional
            The unit cell's box that crystal vectors are taken with respect to.
            If not given, will use a cubic box with a=1 meaning that burgers,
            ξ_uvw and slip_hkl will be interpreted as Cartesian vectors.
        m : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's m-axis,
            i.e. the in-plane direction perpendicular to the dislocation line.  Also
            accepts str values of 'x', 'y', or 'z', in which case the dislocation axis will
            be aligned with the corresponding Cartesian axis.  Default value is 'x'.
        n : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's n-axis,
            i.e. the slip plane normal. Also accepts str values of 'x', 'y', or 'z', in
            which case the dislocation axis will be aligned with the corresponding Cartesian
            axis. Default value is 'y'.
        cart_axes : bool, optional
            Setting this to True will also perform an assertion check that the m- and n-axes
            are both aligned with Cartesian axes. This is a requirement for some of the
            atomic configuration generators. Default value is False as the elastic solution
            by itself does not require the limitation.
        tol : float, optional
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        self.solve(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl, transform=transform,
                   axes=axes, box=box, m=m, n=n, cart_axes=cart_axes, tol=tol)

    def solve(self,
              C: ElasticConstants,
              burgers: npt.ArrayLike,
              ξ_uvw: Optional[npt.ArrayLike] = None,
              slip_hkl: Optional[npt.ArrayLike] = None,
              transform: Optional[npt.ArrayLike] = None,
              axes: Optional[npt.ArrayLike] = None,
              box: Optional[Box] = None,
              m: Union[str, npt.ArrayLike] = 'x',
              n: Union[str, npt.ArrayLike] = 'y',
              cart_axes: bool = False,
              tol: float = 1e-8):
        """
        Computes parameters for a dislocation solution.
        !!!!GENERIC CLASS ONLY HANDLES INPUT PARAMETERS AND DOES NOT SOLVE!!!!

        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Burgers vector.
        ξ_uvw : array-like object
            The Miller crystal vector associated with the dislocation's line
            direction.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        slip_hkl : array-like object
            The Miller plane indices associated with the dislocation's slip
            plane.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        transform : array-like object, optional
            A 3x3 set of orthogonal Cartesian vectors that define the
            transformation matrix to use on C and burgers to convert from the
            standard (unit cell) and dislocation orientations.  The 3 vectors
            will automatically be converted into unit vectors.  Using this is
            an alternative to using ξ_uvw and slip_hkl.
        axes : array-like object, optional
            Same as transform.  Retained for backwards compatibility.
        box : atomman.Box, optional
            The unit cell's box that crystal vectors are taken with respect to.
            If not given, will use a cubic box with a=1 meaning that burgers,
            ξ_uvw and slip_hkl will be interpreted as Cartesian vectors.
        m : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's m-axis,
            i.e. the in-plane direction perpendicular to the dislocation line.  Also
            accepts str values of 'x', 'y', or 'z', in which case the dislocation axis will
            be aligned with the corresponding Cartesian axis.  Default value is 'x'.
        n : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's n-axis,
            i.e. the slip plane normal. Also accepts str values of 'x', 'y', or 'z', in
            which case the dislocation axis will be aligned with the corresponding Cartesian
            axis. Default value is 'y'.
        cart_axes : bool, optional
            Setting this to True will also perform an assertion check that the m- and n-axes
            are both aligned with Cartesian axes. This is a requirement for some of the
            atomic configuration generators. Default value is False as the elastic solution
            by itself does not require the limitation.
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        # Convert burgers, m, n to numpy arrays if needed
        burgers = np.asarray(burgers, dtype=float)

        # Set default box
        if box is None:
            box = Box()

        # Check m, n values
        m, n = self.__mn_check(m, n, cart_axes, tol)

        # ξ_uvw-based orientation parameters
        if ξ_uvw is not None or slip_hkl is not None:
            assert ξ_uvw is not None, 'ξ_uvw and slip_hkl must be given together'
            assert slip_hkl is not None, 'ξ_uvw and slip_hkl must be given together'
            assert transform is None, 'transform cannot be given with ξ_uvw, slip_hkl'
            assert axes is None, 'axes cannot be given with ξ_uvw, slip_hkl'

            transform = self.__find_transform(ξ_uvw, slip_hkl, m, n, box)

        # transform-based orientation parameters
        else:
            if axes is not None:
                assert transform is None, 'axes and transform cannot both be given'
                transform = axes
            if transform is not None:
                transform = axes_check(transform)
            else:
                transform = np.eye(3, dtype=float)

        # Transform burgers and C
        burgers = miller.vector_crystal_to_cartesian(burgers, box)
        burgers = transform.dot(burgers)
        burgers[np.isclose(burgers / np.abs(burgers).max(), 0.0, atol = tol)] = 0.0
        C = C.transform(transform)

        self.__C = C
        self.__m = m
        self.__n = n
        self.__ξ = np.cross(m, n)
        self.__burgers = burgers
        self.__tol = tol
        self.__transform = transform

    def __mn_check(self,
                   m: Union[str, npt.ArrayLike],
                   n: Union[str, npt.ArrayLike],
                   cart_axes: bool,
                   tol: float):
        """Check and verify values for the dislocation solution m- and n-axes."""

        def axis_value(axis):
            """Interpret and validate an axis value"""

            if isinstance(axis, str):

                # Replace str x, y, z values with the corresponding Cartesian vector
                if axis == 'x':
                    axis = np.array([1.0, 0.0, 0.0])
                elif axis == 'y':
                    axis = np.array([0.0, 1.0, 0.0])
                elif axis == 'z':
                    axis = np.array([0.0, 0.0, 1.0])
            else:

                # Convert axis to numpy array and validate
                axis = np.asarray(axis, dtype=float)
                assert axis.shape == (3, ), "m, n axes must be 3D vectors"
                assert np.isclose(np.linalg.norm(m), 1.0, atol=tol, rtol=0.0), "m, n axes must be unit vectors"

                # Check if axis is aligned with a Cartesian axis
                if cart_axes:
                    assert np.isclose(axis, 1.0, atol=tol, rtol=0.0).sum() == 1, "m, n axes must be aligned with Cartesian axes"

            return axis

        # Check m, n independently
        m = axis_value(m)
        n = axis_value(n)

        # Check that they are perpendicular
        assert np.isclose(np.dot(m, n), 0.0, atol=tol, rtol=0.0), "m, n axes must be perpendicular to each other"

        return m, n

    def __find_transform(self, 
                         ξ_uvw: npt.ArrayLike,
                         slip_hkl: npt.ArrayLike,
                         m: npt.ArrayLike,
                         n: npt.ArrayLike,
                         box: Box) -> np.ndarray:
        """
        Find the transformation matrix to the dislocation system based on
        the crystallographic slip plane and line vector.
        """

        # ξ_axis is the Cartesian unit vector along the dislocation line
        ξ_axis = box.vector_crystal_to_cartesian(ξ_uvw)
        ξ_axis = ξ_axis / np.linalg.norm(ξ_axis)

        # n_axis is the unit vector normal to slip plane
        n_axis = box.plane_crystal_to_cartesian(slip_hkl)

        # m_axis is the in-plane unit vector perpendicular to ξ
        m_axis = np.cross(n_axis, ξ_axis)

        # Construct the transformation matrix
        transform = np.array([m_axis, n_axis, ξ_axis])

        # Transform transform to m, n, ξ orientation
        T = np.array([m, n, np.cross(m, n)]).T

        transform = T.dot(transform)

        return transform

    @property
    def m(self) -> np.ndarray:
        """numpy.ndarray: The Cartesian vector for orienting the dislocation system's m vector"""
        return self.__m

    @property
    def n(self) -> np.ndarray:
        """numpy.ndarray: The Cartesian vector for orienting the dislocation system's n vector"""
        return self.__n

    @property
    def ξ(self) -> np.ndarray:
        """numpy.ndarray: The Cartesian vector for orienting the dislocation system's line vector"""
        return self.__ξ

    @property
    def C(self) -> ElasticConstants:
        """atomman.ElasticConstants: The associated elastic constants"""
        return self.__C

    @property
    def burgers(self) -> np.ndarray:
        """numpy.ndarray: The Cartesian Burgers vector"""
        return self.__burgers

    @property
    def tol(self) -> float:
        """float: The tolerance value used to verify the solution"""
        return self.__tol

    @property
    def transform(self) -> np.ndarray:
        """numpy.ndarray: The transformation matrix associated with the dislocation system's orientation"""
        return self.__transform

    def characterangle(self, unit: str = 'degree') -> float:
        """
        Returns the dislocation's character angle.

        Parameters
        ----------
        unit : str, optional
            Specify whether the angle is given in 'degree' (default)
            or 'radian'.

        Returns
        -------
        float
            The dislocation character angle.
        """
        return vect_angle(self.burgers, self.ξ, unit=unit)

    @property
    def K_coeff(self) -> float:
        """float : The energy coefficient"""

        # K = b_i K_ij b_j / (b_k b_k)
        return (self.burgers.dot(self.K_tensor.dot(self.burgers))
                / self.burgers.dot(self.burgers))

    @property
    def K_tensor(self) -> np.ndarray:
        """numpy.ndarray : The energy coefficient tensor"""
        raise NotImplementedError('Needs to be defined by subclass')

    @property
    def preln(self) -> float:
        """float : The pre-ln strain energy factor"""
        # a = b_i K_ij b_j / (4 π)
        return self.burgers.dot(self.K_tensor.dot(self.burgers)) / (4 * np.pi)

    def displacement(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent displacements.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3D vector displacements at all given points.
        """
        raise NotImplementedError('Needs to be defined by subclass')

    def strain(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent strains.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 strain states at all given points.
        """
        raise NotImplementedError('Needs to be defined by subclass')

    def stress(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent stresses.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 stress states at all given points.
        """
        raise NotImplementedError('Needs to be defined by subclass')
