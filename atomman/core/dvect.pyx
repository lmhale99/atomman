# coding: utf-8
# cython: language_level=3
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# http://cython.org/
import cython

# http://www.numpy.org/
import numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
def dvect(pos_0, pos_1, box, pbc):
    """
    Computes the shortest vector between pos_0 and pos_1 using box
    dimensions and accounting for periodic boundaries.

    Parameters
    ----------
    pos_0 : numpy.ndarray, list, or tuple
        Absolute Cartesian vector position(s) to use as reference point(s).
    pos_1 : numpy.ndarray, list, or tuple
        Absolute Cartesian vector position(s) to find relative to pos_0.
    box : atomman.Box
        Defines the system/box dimensions
    pbc : list, tuple, or numpy.ndarray of bool.
        Three Boolean values indicating which of the three box vectors are
        periodic (True means periodic).

    Returns
    -------
    numpy.ndarray
        The shortest vectors from each pos_0 to pos_1 positions.
    """
    
    # Convert pos_0 to numpy array with proper dimensions
    pos_0 = np.asarray(pos_0, dtype=np.float64)
    if pos_0.ndim == 0:
        raise TypeError('Invalid pos_0')
    if pos_0.ndim == 1:
        pos_0 = pos_0[np.newaxis, :]

    # Convert pos_1 to numpy array with proper dimensions
    pos_1 = np.asarray(pos_1, dtype=np.float64)
    if pos_1.ndim == 0:
        raise TypeError('Invalid pos_1')
    if pos_1.ndim == 1:
        pos_1 = pos_1[np.newaxis, :]

    # Broadcast to compatible lengths
    if len(pos_0) == 1:
        pos_0 = np.broadcast_to(pos_0, pos_1.shape)
    elif len(pos_1) == 1:
        pos_1 = np.broadcast_to(pos_1, pos_0.shape)
    elif len(pos_0) != len(pos_1):
        raise ValueError('Incompatible pos lengths')

    # Extract box vectors
    bvects = box.vects

    # Call the cython function
    return dvect_c(pos_0, pos_1, bvects, pbc[0], pbc[1], pbc[2])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef dvect_c(const double[:,:] pos_0,
             const double[:,:] pos_1, 
             const double[:,:] bvects,
             const bint pbc_x,
             const bint pbc_y,
             const bint pbc_z):
    """
    Computes the shortest distance between pos_0 and pos_1 using box 
    dimensions and accounting for periodic boundaries.
    
    Parameters
    ----------
    pos_0 : cython.memoryview
        Absolute Cartesian vector position(s) to use as reference point(s).
    pos_1 : cython.memoryview
        Absolute Cartesian vector position(s) to find relative to pos_0.
    bvects : cython.memoryview
        3x3 array defining the system/box dimensions.
    pbc_x : bint
        Flag indicating to make x periodic.
    pbc_y : bint
        Flag indicating to make y periodic.
    pbc_z : bint
        Flag indicating to make z periodic.

    Returns
    -------
    cython.memoryview
        The shortest vectors from each pos_0 to pos_1 positions.
    """
    
    # Define parameters
    cdef Py_ssize_t ni = pos_0.shape[0]
    cdef int nj = 3
    cdef int i, j, x, y, z, xl, xh, yl, yh, zl, zh    
    cdef double[:] test = np.empty(3)
    cdef double mag_test, mag_d
    
    # Define output array and its view
    d = np.empty_like(pos_0, np.float64)
    cdef double[:,:] dv = d
    
    # Create iterators based on pbc
    if pbc_x:
        xl, xh = -1, 2
    else:
        xl, xh = 0, 1
    if pbc_y:
        yl, yh = -1, 2
    else:
        yl, yh = 0, 1
    if pbc_z:
        zl, zh = -1, 2
    else:
        zl, zh = 0, 1
    
    # Loop over all pos
    for i in range(ni):
       
        # Compute pos_1 - pos_0 
        for j in range(nj):
            dv[i,j] = pos_1[i,j] - pos_0[i,j]
        
        # Loop over all periodic boundary conditions
        for x in range(xl, xh):
            for y in range(yl, yh):
                for z in range(zl, zh):
                    if x == 0 and y == 0 and z == 0:
                        continue
                    
                    # Compute pos_1 - pos_0 + boundary image shifts
                    for j in range(nj):
                        test[j] = (pos_1[i,j] - pos_0[i,j] 
                                   + x * bvects[0,j] 
                                   + y * bvects[1,j] 
                                   + z * bvects[2,j])
                    
                    # Replace d if new vector is smaller
                    mag_test = test[0] * test[0] + test[1] * test[1] + test[2] * test[2]
                    mag_d = dv[i,0] * dv[i,0] + dv[i,1] * dv[i,1] + dv[i,2] * dv[i,2]
                    if mag_test < mag_d:
                        for j in range(nj):
                            dv[i,j] = test[j]

    return d