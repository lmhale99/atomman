# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
                        
# http://cython.org/
import cython
cimport cython

# http://www.numpy.org/
import numpy as np
cimport numpy as np

def dvect_cython(pos_0, pos_1, box, pbc):
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
    """
    # Convert pos_0 to numpy array with proper dimensions
    pos_0 = np.asarray(pos_0)
    if pos_0.ndim == 1:
        pos_0 = pos_0[np.newaxis, :]
    
    # Convert pos_1 to numpy array with proper dimensions
    pos_1 = np.asarray(pos_1)
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
    return cdvect(pos_0, pos_1, bvects, pbc)

@cython.boundscheck(False)
cdef np.ndarray[np.float64_t, ndim=2] cdvect(
                                             np.ndarray[np.float64_t, ndim=2] pos_0, 
                                             np.ndarray[np.float64_t, ndim=2] pos_1, 
                                             np.ndarray[np.float64_t, ndim=2] bvects, 
                                             pbc):
    """
    Computes the shortest distance between pos_0 and pos_1 using box 
    dimensions and accounting for periodic boundaries.
    
    Parameters
    ----------
    pos_0 : numpy.ndarray
        Absolute Cartesian vector position(s) to use as reference point(s).
    pos_1 : numpy.ndarray
        Absolute Cartesian vector position(s) to find relative to pos_0.
    bvects : numpy.ndarray
        3x3 array defining the system/box dimensions.
    pbc : list, tuple, or numpy.ndarray of bool.
        Three Boolean values indicating which of the three box vectors are
        periodic (True means periodic).
    """
    
    # Define parameters
    cdef Py_ssize_t ni = pos_0.shape[0]
    cdef Py_ssize_t nj = 3
    
    cdef Py_ssize_t i, j, x, y, z, xl, xu, yl, yu, zl, zu
        
    cdef np.ndarray[np.float64_t, ndim=2] d = np.empty_like(pos_0)
    cdef np.ndarray[np.float64_t] test = np.empty(3)
    cdef np.float64_t test2, d2
    
    # Create iterators based on pbc
    if pbc[0] is True:
        xl, xh = -1, 2
    else:
        xl, xh = 0, 1
    if pbc[1] is True:
        yl, yh = -1, 2
    else:
        yl, yh = 0, 1
    if pbc[2] is True:
        zl, zh = -1, 2
    else:
        zl, zh = 0, 1
    
    # Loop over all pos
    for i in range(ni):
        # Compute pos_1 - pos_0 
        for j in range(nj):
            d[i,j] = pos_1[i,j] - pos_0[i,j] 
        
        # Loop over all periodic boundary conditions
        for x in range(xl, xh):
            for y in range(yl, yh):
                for z in range(zl, zh):
                    
                    # Compute pos_1 - pos_0 + boundary image shifts
                    for j in range(nj):
                        test[j] = (pos_1[i,j] - pos_0[i,j] 
                                   + x * bvects[0,j] 
                                   + y * bvects[1,j] 
                                   + z * bvects[2,j])
                    
                    # Replace d if new vector is smaller
                    test2 = test[0]*test[0]+test[1]*test[1]+test[2]*test[2]
                    d2 = d[i,0]*d[i,0]+d[i,1]*d[i,1]+d[i,2]*d[i,2]
                    if test2 < d2:
                        for j in range(nj):
                            d[i,j] = test[j]
    
    return d