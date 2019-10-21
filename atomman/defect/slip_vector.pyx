# coding: utf-8
# cython: language_level=3
# Standard Python libraries

# http://cython.org/
import cython

# http://www.numpy.org/
import numpy as np

# atomman imports
from .. import NeighborList
from ..core.dvect cimport dvect_c

@cython.boundscheck(False)
@cython.wraparound(False)
def slip_vector(system_0, system_1, neighbors=None, cutoff=None):
    """
    Compute the slip vectors for all atoms.  Note that this differs from the
    original formulation in that it is not normalized by number of slipped
    neighbors.
    
        s_i = -Σ_j d_ij(t) - d_ij(0)
    
    where j is neighbor atoms of atom i, and d_ij() is vector distance between
    atoms i and j at time t.
    
    Parameters
    ----------
    system_0 : atomman.system
        The base/reference system to use.
    system_1 : atomman.system
        The defect/current system to use.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system_0 to use.  Either neighbors
        or cutoff must be given, or system_0 must have a neighbors attribute.
    cutoff : float, optional
        Cutoff distance for computing a neighbor list for system_0.  Either
        neighbors or cutoff must be given, or system_0 have a neighbors
        attribute.
    
    Returns
    -------
    numpy.ndarray
        The computed slip vectors.
    """
    # Check that the two systems have the same number of atoms
    if system_0.natoms != system_1.natoms:
        raise ValueError('systems have different number of atoms')

    # Neighbor list setup
    if neighbors is not None:
        assert cutoff is None, 'neighbors and cutoff cannot both be given'
    elif cutoff is not None:
        neighbors = NeighborList(system=system_0, cutoff=cutoff)
    elif hasattr(system_0, 'neighbors'):
        neighbors = system_0.neighbors
    else:
        raise ValueError('neighbors or cutoff is required')
    
    pos_0 = system_0.atoms.pos
    pos_1 = system_1.atoms.pos
    bvects = system_0.box.vects
    nlist = neighbors.nlist
    
    return slip_vector_c(pos_0, pos_1, bvects, nlist, system_0.pbc[0], system_0.pbc[1], system_0.pbc[2])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef slip_vector_c(const double[:,:] pos_0, 
                   const double[:,:] pos_1, 
                   const double[:,:] bvects, 
                   const long long [:,:] nlist, 
                   bint pbc_a, bint pbc_b, bint pbc_c):
    
    cdef Py_ssize_t i, j, coordmax, ni, n, coord
    
    cdef double[:,:] d_0, d_1
    
    # Define the slip vector
    slip = np.zeros((pos_0.shape[0], 3))
    cdef double[:,:] slipv = slip
    
    coordmax = 0
    for i in range(nlist.shape[0]):
        if nlist[i, 0] > coordmax:
            coordmax = nlist[i, 0]
    
    cdef double[:,:] upos_0, vpos_0, upos_1, vpos_1
    upos_0 = np.zeros((coordmax, 3))
    vpos_0 = np.zeros((coordmax, 3))
    upos_1 = np.zeros((coordmax, 3))
    vpos_1 = np.zeros((coordmax, 3))
    
    for i in range(pos_0.shape[0]):
        coord = nlist[i, 0]
        for n in range(coord):
            for j in range(3):
                ni = nlist[i, n+1]
                upos_0[n, j] = pos_0[i, j]
                vpos_0[n, j] = pos_0[ni, j]
                upos_1[n, j] = pos_1[i, j]
                vpos_1[n, j] = pos_1[ni, j]
                
        d_0 = dvect_c(upos_0, vpos_0, bvects, pbc_a, pbc_b, pbc_c)
        d_1 = dvect_c(upos_1, vpos_1, bvects, pbc_a, pbc_b, pbc_c)
        
        for n in range(coord):
            for j in range(3):
                slipv[i, j] -= d_1[n, j] - d_0[n, j]
    
    return slip