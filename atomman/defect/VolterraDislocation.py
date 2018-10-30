# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..tools import axes_check, vect_angle

class VolterraDislocation(object):
    """
    Generic class for a Volterra solution of a straight dislocation.
    """
    
    def __init__(self, C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
        """
        Initializes the solution. Calls solve.
        
        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Cartesian Burgers vector.
        axes : array-like object, optional
            3x3 set of rotational axes for the system. If given, burgers
            will be transformed using axes.
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
        self.solve(C, burgers, axes=axes, m=m, n=n, tol=tol)
        
    def solve(self, C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
        """
        Computes parameters for a dislocation solution.
        !!!!NOT IMPLEMENTED FOR GENERIC CLASS!!!!
        
        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Cartesian Burgers vector.
        axes : array-like object, optional
            3x3 set of rotational axes for the system. If given, burgers
            will be transformed using axes.
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
        burgers = np.asarray(burgers, dtype='float64')
        m = np.asarray(m, dtype='float64')
        n = np.asarray(n, dtype='float64')
        
        # Check m, n values
        assert np.isclose(np.linalg.norm(m), 1.0, atol = tol), "m must be a unit vector"
        assert np.isclose(np.linalg.norm(n), 1.0, atol = tol), "n must be a unit vector"
        assert np.isclose(np.dot(m, n), 0.0, atol = tol), "m and n must be perpendicular"
        
        # Transform burgers and C
        if axes is not None:
            T = axes_check(axes)
            burgers = T.dot(burgers)
            burgers[np.isclose(burgers / burgers.max(), 0.0, atol = tol)] = 0.0
            C = C.transform(axes)
        
        self.__C = C
        self.__m = m
        self.__n = n
        self.__u = np.cross(m, n)
        self.__burgers = burgers
        self.__tol = tol
    
    @property
    def m(self):
        return self.__m
    
    @property
    def n(self):
        return self.__n
    
    @property
    def u(self):
        return self.__u
    
    @property
    def C(self):
        return self.__C
    
    @property
    def burgers(self):
        return self.__burgers
    
    @property
    def tol(self):
        return self.__tol
    
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
        return vect_angle(self.__burgers, self.__u, unit=unit)
    
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