# coding: utf-8
# Standard Python libraries
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from . import dislocation_system_transform
from .. import Box
from ..tools import axes_check, vect_angle, miller

class VolterraDislocation(object):
    """
    Generic class for a Volterra solution of a straight dislocation.
    """
    
    def __init__(self, C, burgers, ξ_uvw=None, slip_hkl=None, transform=None,
                 axes=None, box=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
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
            standard (unit cell) and dislocation orientations.  The 3 vectors
            will automatically be converted into unit vectors.  Using this is
            an alternative to using ξ_uvw and slip_hkl.
        axes : array-like object, optional
            Same as transform.  Retained for backwards compatibility.
        box : atomman.Box, optional
            The unit cell's box that crystal vectors are taken with respect to.
            If not given, will use a cubic box with a=1 meaning that burgers,
            ξ_uvw and slip_hkl will be interpreted as Cartesian vectors.
        m : array-like object, optional
            The m unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [1,0,0]
            (x-axis).
        n : array-like object, optional
            The n unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [0,1,0]
            (y-axis). n is normal to the dislocation slip plane.
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        self.solve(C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl, transform=transform,
                   axes=axes, box=box, m=m, n=n, tol=tol)
        
    def solve(self, C, burgers, ξ_uvw=None, slip_hkl=None, transform=None,
              axes=None, box=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
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
        m : array-like object, optional
            The m unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [1,0,0]
            (x-axis).
        n : array-like object, optional
            The n unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [0,1,0]
            (y-axis). n is normal to the dislocation slip plane.
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        # Convert burgers, m, n to numpy arrays if needed
        burgers = np.asarray(burgers, dtype=float)
        m = np.asarray(m, dtype=float)
        n = np.asarray(n, dtype=float)
        
        # Check m, n values
        assert np.isclose(np.linalg.norm(m), 1.0, atol = tol), "m must be a unit vector"
        assert np.isclose(np.linalg.norm(n), 1.0, atol = tol), "n must be a unit vector"
        assert np.isclose(np.dot(m, n), 0.0, atol = tol), "m and n must be perpendicular"
        
        # ξ_uvw-based orientation parameters
        if ξ_uvw is not None or slip_hkl is not None:
            assert ξ_uvw is not None, 'ξ_uvw and slip_hkl must be given together'
            assert slip_hkl is not None, 'ξ_uvw and slip_hkl must be given together'
            assert transform is None, 'transform cannot be given with ξ_uvw, slip_hkl'
            assert axes is None, 'axes cannot be given with ξ_uvw, slip_hkl'
            
            transform = dislocation_system_transform(ξ_uvw, slip_hkl, m=m, n=n, box=box, tol=tol)
        
        # transform-based orientation parameters
        else:
            if axes is not None:
                assert transform is None, 'axes and transform cannot both be given'
                transform = axes
            if transform is not None:
                transform = axes_check(transform)
            else:
                transform = np.eye(3, dtype=float)

        # Set default box
        if box is None:
            box = Box()
        
        # Transform burgers and C
        burgers = miller.vector_crystal_to_cartesian(burgers, box)
        burgers = transform.dot(burgers)
        burgers[np.isclose(burgers / burgers.max(), 0.0, atol = tol)] = 0.0
        C = C.transform(transform)
        
        self.__C = C
        self.__m = m
        self.__n = n
        self.__ξ = np.cross(m, n)
        self.__burgers = burgers
        self.__tol = tol
        self.__transform = transform
    
    @property
    def m(self):
        return self.__m
    
    @property
    def n(self):
        return self.__n
    
    @property
    def ξ(self):
        return self.__ξ
    
    @property
    def C(self):
        return self.__C
    
    @property
    def burgers(self):
        return self.__burgers
    
    @property
    def tol(self):
        return self.__tol
    
    @property
    def transform(self):
        return self.__transform

    def characterangle(self, unit='degree'):
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
    def K_coeff(self):
        """float : The energy coefficient"""
        
        # K = b_i K_ij b_j / (b_k b_k)
        return (self.burgers.dot(self.K_tensor.dot(self.burgers))
                / self.burgers.dot(self.burgers))
    
    @property
    def K_tensor(self):
        raise NotImplementedError('Needs to be defined by subclass')
    
    @property
    def preln(self):
        """float : The pre-ln strain energy factor"""
        # a = b_i K_ij b_j / (4 π)
        return self.burgers.dot(self.K_tensor.dot(self.burgers)) / (4 * np.pi)
    
    def displacement(self, pos):
        raise NotImplementedError('Needs to be defined by subclass')
    
    def stress(self, pos):
        raise NotImplementedError('Needs to be defined by subclass')