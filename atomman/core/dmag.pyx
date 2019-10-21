# coding: utf-8
# cython: language_level=3
# Standard Python libraries
from copy import deepcopy

# http://cython.org/
import cython

# http://www.numpy.org/
import numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
def dmag(pos_0, pos_1, box, pbc):
    """
    Computes the shortest distance between pos_0 and pos_1 using box
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
        The shortest vector magnitude from each pos_0 to pos_1 positions.
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
    return dmag2_c(pos_0, pos_1, bvects, pbc[0], pbc[1], pbc[2])**0.5

@cython.boundscheck(False)
@cython.wraparound(False)
cdef dmag2_c(const double[:,:] pos_0,
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
        The shortest vector magnitude squared from each pos_0 to pos_1 positions.
    """
    
    # Define parameters
    cdef Py_ssize_t ni = pos_0.shape[0]
    cdef Py_ssize_t nj = 3
    cdef Py_ssize_t i, j, x, y, z, xl, xh, yl, yh, zl, zh
    cdef double mag2_test
    cdef double[:] d = np.empty(3, dtype=np.float64)
    
    # Define output array and its view
    mag2_d = np.empty(ni, dtype=np.float64)
    cdef double [:] mag2_dv = mag2_d

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
            d[j] = pos_1[i,j] - pos_0[i,j]
        mag2_dv[i] = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        
        # Loop over all periodic boundary conditions
        for x in range(xl, xh):
            for y in range(yl, yh):
                for z in range(zl, zh):
                    if x == 0 and y == 0 and z == 0:
                        continue
                    
                    # Compute pos_1 - pos_0 + boundary image shifts
                    for j in range(nj):
                        d[j] = (pos_1[i,j] - pos_0[i,j] 
                                + x * bvects[0,j] 
                                + y * bvects[1,j] 
                                + z * bvects[2,j])
                    
                    # Replace d if new vector is smaller
                    mag2_test = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
                    if mag2_test < mag2_dv[i]:
                        mag2_dv[i] = mag2_test

    return mag2_d 