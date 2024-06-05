# coding: utf-8
# cython: language_level=3

# http://cython.org/
import cython

# http://www.numpy.org/
import numpy as np

# atomman imports
from .dmag cimport dmag2_c   

@cython.boundscheck(False)
@cython.wraparound(False)
def nlist(system, double cutoff, Py_ssize_t initialsize=20, Py_ssize_t deltasize=10):
    """
    Calculates a neighbor list for all atoms in a System taking periodic
    boundaries into account.
    
    Parameters
    ----------
    system : atomman.System 
        The system to calculate the neighbor list for.
    cutoff : float
        Radial cutoff distance for identifying neighbors.
    initialsize : int, optional
        The number of neighbor positions to initially assign to each atom.
        Default value is 20.
    deltasize : int, optional
        Specifies the number of extra neighbor positions to allow each atom
        when the number of neighbors exceeds the underlying array size.
        Default value is 10.

    Returns
    -------
    numpy.ndarray of int
        Array listing number of neighbors and neighbor ids for each atom in
        System.  First term in each row is the atom's coordination number, c.
        The next c values are the atom's neighbor ids.  
    """
    
    # Define variables based on input parameters
    pos = system.atoms.pos
    cdef const double[:,:] posv = pos
    cdef const double[:,:] vects = system.box.vects
    cdef const double[:] origin = system.box.origin
    cdef bint pbc_a = system.pbc[0]
    cdef bint pbc_b = system.pbc[1]
    cdef bint pbc_c = system.pbc[2]
    cdef Py_ssize_t maxneighbors = initialsize
    cdef Py_ssize_t maxatomsperbin = 40
    cdef Py_ssize_t natoms = posv.shape[0]
    cdef double cutoff2 = cutoff*cutoff
   
    # Define basic iteration index variables
    cdef Py_ssize_t i, j, k, l

    # Define superbox identification parameters
    cdef Py_ssize_t x, y, z
    cdef double corner
    cdef double[:] supermin=np.empty(3)
    cdef double[:] supermax=np.empty(3)
    
    # Define bins
    cdef double binsize
    cdef double[:] xbins, ybins, zbins
    cdef Py_ssize_t numxbins, numybins, numzbins
    
    # Define bin indexers and associated atom indexes
    #cdef long long [:] xindex, yindex, zindex
    #cdef long long [:,:] xyzindex
    cdef long long [:] atomindex
    cdef long long [:,:] newxyzindex
    cdef long long [:] newatomindex
    
    # Define ghost atom info
    cdef Py_ssize_t xl, xh, yl, yh, zl, zh
    cdef double[:,:] ghostpos = np.empty((0, 3))
    cdef double[:,:] newghostpos
    newpos = np.empty(pos.shape)
    cdef double[:,:] newposv = newpos
    cdef long long[:] ghostindex = np.empty(0, dtype=np.int64)
    cdef long long[:] newindex = np.empty(posv.shape[0], dtype=np.int64)
    cdef long long[:] newghostindex 
    cdef long long [:,:] xyzghostindex
    
    # Define final bins
    cdef Py_ssize_t maxc, c, n
    cdef Py_ssize_t dc, dx, dy, dz
    cdef long long[:, :, :, :] xyzbins, newbins
    
    # Define lists of atoms to compare
    cdef Py_ssize_t uindex, vindex, u, v, w
    cdef long long [:] shortlist, longlist, superlonglist
    cdef bint end
    
    # Define positions and distances between them
    cdef double[:,:] upos, vpos
    cdef double[:] dmag2
    
    # Define neighbors list and assignment terms
    cdef Py_ssize_t uj, vj
    cdef long long[:, :] neighbors = np.empty((natoms, maxneighbors+1), dtype=np.int64)
    for i in range(natoms):
        neighbors[i, 0] = 0
    cdef long long[:, :] newneighbors
    cdef bint isnew

    # Determine orthogonal superbox that fully encompasses the system
    for j in range(3):
        supermin[j] = origin[j]
        supermax[j] = origin[j]
    for z in range(0, 2):
        for y in range(0, 2):
            for x in range(0, 2):
                for j in range(3):
                    corner = origin[j] + x * vects[0, j] + y * vects[1, j] + z * vects[2, j]
                    if corner < supermin[j]:
                        supermin[j] = corner
                    if corner > supermax[j]:
                        supermax[j] = corner
    for j in range(3):
        supermin[j] -= 1.01 * cutoff
        supermax[j] += 1.01 * cutoff
    
    # Construct bins
    binsize = cutoff
    xbins = np.arange(supermin[0], supermax[0] + binsize, binsize)
    ybins = np.arange(supermin[1], supermax[1] + binsize, binsize)
    zbins = np.arange(supermin[2], supermax[2] + binsize, binsize)
    numxbins = len(xbins)
    numybins = len(ybins)
    numzbins = len(zbins)
    
    # Build xyz box index for each atom
    xindex = np.digitize(pos[:, 0], xbins) - 1
    yindex = np.digitize(pos[:, 1], ybins) - 1
    zindex = np.digitize(pos[:, 2], zbins) - 1
    xyzindex = np.hstack((xindex[:, np.newaxis], yindex[:, np.newaxis], zindex[:, np.newaxis]))
    
    # Relate atom's id to xyz_index
    atomindex = np.arange(natoms, dtype=np.int64)
    
    # Identify all bins with real atoms
    realbins = unique_rows2(xyzindex)
    
    # Create iterators based on pbc
    if pbc_a:
        xl, xh = -1, 2
    else:
        xl, xh = 0, 1
    if pbc_b:
        yl, yh = -1, 2
    else:
        yl, yh = 0, 1
    if pbc_c:
        zl, zh = -1, 2
    else:
        zl, zh = 0, 1

    # Construct list of ghost atoms in the superbox
    for x in range(xl, xh):
        for y in range(yl, yh):
            for z in range(zl, zh):
                if x == 0 and y == 0 and z == 0:
                    pass
                else:
                    k=0
                    for i in range(posv.shape[0]):
                        for j in range(3):
                            newposv[i, j] = x * vects[0, j] + y * vects[1, j] + z * vects[2, j] + posv[i, j]

                        if (    newposv[i, 0] > supermin[0] and newposv[i, 0] < supermax[0]
                            and newposv[i, 1] > supermin[1] and newposv[i, 1] < supermax[1]
                            and newposv[i, 2] > supermin[2] and newposv[i, 2] < supermax[2]):
                            
                            newindex[k] = i
                            k += 1

                    ghostpos = np.vstack((ghostpos, newpos[newindex[:k]]))
                    ghostindex = np.hstack((ghostindex, newindex[:k]))

    # Append xyzindex and atomindex lists with ghost atoms
    if len(ghostpos) > 0:
        xindex = np.digitize(ghostpos[:, 0], xbins) - 1
        yindex = np.digitize(ghostpos[:, 1], ybins) - 1
        zindex = np.digitize(ghostpos[:, 2], zbins) - 1
        xyzghostindex = np.hstack((xindex[:, np.newaxis], 
                                   yindex[:, np.newaxis],
                                   zindex[:, np.newaxis]))
        xyzindex = np.vstack((xyzindex, xyzghostindex))
        atomindex = np.hstack((atomindex, ghostindex))
    
    # Assign atoms and ghost atoms to xyz bins
    maxc = 0
    xyzbins = np.zeros((numxbins, numybins, numzbins, maxatomsperbin + 1), dtype=np.int64)
    for n in range(atomindex.shape[0]):
        x, y, z = xyzindex[n]
        c = xyzbins[x, y, z, 0] + 1
        
        # Increase size if needed
        if c == maxatomsperbin:
            newbins = np.zeros((numxbins, numybins, numzbins, maxatomsperbin + 11), dtype=np.int64)
            for i in range(xyzbins.shape[0]):
                for j in range(xyzbins.shape[1]):
                    for k in range(xyzbins.shape[2]):
                        for l in range(maxatomsperbin + 1):
                            newbins[i, j, k, l] = xyzbins[i, j, k, l]
            xyzbins = newbins
            maxatomsperbin += 10
        
        if c > maxc:
            maxc = c
        xyzbins[x, y, z, 0] = c
        xyzbins[x, y, z, c] = atomindex[n]

    superlonglist = np.empty(14 * maxc, dtype=np.int64)

    # Iterate over all bins with real atoms
    for i in range(len(realbins)):
        x, y, z = realbins[i]
        c = xyzbins[x, y, z, 0]
        shortlist = np.empty(c, dtype=np.int64)
        
        # shortlist = all atoms in current bin
        for j in range(c):
            shortlist[j] = xyzbins[x, y, z, j+1]
            superlonglist[j] = shortlist[j]
        
        # Add all atoms in half of the nearby bins to longlist
        end = False
        for dz in range(-1, 2):
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    
                    # Stop when reach center bin
                    if dx == 0 and dy == 0 and dz == 0:
                        end = True
                        break
                    
                    # Skip non-existant neighbor bins
                    if (x + dx < 0 or x + dx == numxbins or
                        y + dy < 0 or y + dy == numybins or
                        z + dz < 0 or z + dz == numzbins):
                        continue
                    
                    dc = xyzbins[x + dx, y + dy, z + dz, 0]
                    for j in range(dc):
                        superlonglist[c+j] = xyzbins[x + dx, y + dy, z + dz, j+1]
                    c += dc
                if end:
                    break
            if end:
                break
                
        longlist = superlonglist[:c]
        
        # Compare all atoms in shortlist to longlist.
        for u in range(shortlist.shape[0]):
            uindex = shortlist[u]
            
            upos = np.empty((longlist.shape[0]-u-1, 3))
            vpos = np.empty((longlist.shape[0]-u-1, 3))

            for w, v in enumerate(range(u+1, longlist.shape[0])):
                for j in range(3):
                    vindex = longlist[v]
                    upos[w, j] = posv[uindex, j]
                    vpos[w, j] = posv[vindex, j]
            
            # Compute distances
            dmag2 = dmag2_c(upos, vpos, vects, pbc_a, pbc_b, pbc_c)
            
            # Assign neighbors if within cutoff
            for w, v in enumerate(range(u+1, longlist.shape[0])):
                if dmag2[w] < cutoff2:
                    vindex = longlist[v]
            
                    if uindex != vindex:
                        isnew = True
                        uj = -1
                        vj = -1
                        
                        # Find uj position to insert vindex
                        for j in range(1, neighbors[uindex, 0] + 1):
                            # Check uindex's neighbors for vindex
                            if neighbors[uindex, j] == vindex:
                                isnew = False
                                break
                            elif neighbors[uindex, j] > vindex:
                                uj = j
                                break
                        if uj == -1:
                            uj = neighbors[uindex, 0] + 1

                        if isnew:
                            # Find vj position to insert uindex
                            for j in range(1, neighbors[vindex, 0] + 1):
                                if neighbors[vindex, j] > uindex:
                                    vj = j
                                    break
                            if vj == -1:
                                vj = neighbors[vindex, 0] + 1

                            # Increase coordination
                            neighbors[uindex, 0] += 1
                            neighbors[vindex, 0] += 1

                            # Extend system size if needed
                            if neighbors[uindex, 0] > maxneighbors or neighbors[vindex, 0] > maxneighbors:
                                newneighbors = np.empty((natoms, maxneighbors + deltasize + 1), dtype=np.int64)
                                for j in range(neighbors.shape[0]):
                                    for k in range(maxneighbors + 1):
                                        newneighbors[j, k] = neighbors[j, k]
                                neighbors = newneighbors
                                maxneighbors += deltasize

                            # Shift neighbor indices with higher values
                            for j in range(neighbors[uindex, 0], uj - 1, -1):
                                neighbors[uindex, j] = neighbors[uindex, j - 1]
                            for j in range(neighbors[vindex, 0], vj - 1, -1):
                                neighbors[vindex, j] = neighbors[vindex, j - 1]

                            # Assign neighbors to each other
                            neighbors[uindex, uj] = vindex
                            neighbors[vindex, vj] = uindex

    return np.asarray(neighbors)

@cython.boundscheck(False)
@cython.wraparound(False)
def unique_rows2(a):  
    """Takes two-dimensional array a and returns only the unique rows."""
    return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])