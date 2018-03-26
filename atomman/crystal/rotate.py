# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
import atomman.core.System
import atomman.core.supersize
from .miller import vectortocartesian, vector4to3

def rotate(system, hkls, tol=1e-5):
    """
    Transforms a System representing a periodic crystal cell from a standard
    orientation to a specified orientation. Note: if hexagonal indices are
    given, the vectors will be reduced to the smallest hkl integer
    representation.
    
    Parameters
    ----------
    system : atomman.System
        The system to transform
    hkls : numpy.ndarray of int
        A (3, 3) array of the Miller crystal vectors or a (3, 4) array of
        Miller-Bravais hexagonal crystal vectors to use in transforming the
        system.
    tol : float, optional
        Tolerance parameter used in rounding atomic positions near the
        boundaries to the boundary values.  In box-relative coordinates, any
        atomic positions within tol of 0 or 1 will be rounded to 0 or 1,
        respectively.  Default value is 1e-5.
        
    Returns
    -------
    atomman.System
        A new fully periodic system rotated and transformed according to the
        hkls crystal vectors.
    """
    
    # Check parameters
    try:
        hkls = np.asarray(hkls, dtype='int64')
       
        if hkls.shape == (3, 4):
            hkls = atomman.crystal.miller.vector4to3(hkls)
        assert hkls.shape == (3, 3)
    except:
        raise ValueError('Invalid hkls crystal indices')
    
    # Get natoms and volume of system
    natoms = system.natoms
    volume = np.abs(system.box.avect.dot(np.cross(system.box.bvect,
                                                  system.box.cvect)))
    
    # Convert hkls to Cartesian units and compute new volume and natoms
    newvects = atomman.crystal.miller.vectortocartesian(hkls, box=system.box)
    newvolume = np.abs(newvects[0].dot(np.cross(newvects[1], newvects[2])))
    newnatoms = int(round(newvolume / volume) * natoms)
    
    # Check new values
    if newnatoms == 0:
        raise ValueError('New box has no atoms/volume: vectors are parallel or planar')
    
    # Identify box corners of new system wrt hkls
    corners = np.empty((8,3), dtype='int64')
    corners[0] = np.zeros(3)
    corners[1] = hkls[0]
    corners[2] = hkls[1]
    corners[3] = hkls[2]
    corners[4] = hkls[0] + hkls[1]
    corners[5] = hkls[0] + hkls[2]
    corners[6] = hkls[1] + hkls[2]
    corners[7] = hkls[0] + hkls[1] + hkls[2]
    
    # Create a supercell of system that contains all box corners
    a_mults = (corners[:,0].min()-1, corners[:,0].max()+1)
    b_mults = (corners[:,1].min()-1, corners[:,1].max()+1)
    c_mults = (corners[:,2].min()-1, corners[:,2].max()+1)
    system2 = atomman.core.supersize(system, a_mults, b_mults, c_mults)
    
    # Change system.box.vects to newvects
    system2.box_set(vects=newvects, scale=False)
    
    # Round atom positions near box boundaries to the boundaries
    spos = system2.atoms_prop('pos', scale=True)
    spos[np.isclose(spos, 0.0, atol=tol)] = 0.0
    spos[np.isclose(spos, 1.0, atol=tol)] = 1.0
    
    # Identify all atoms whose positions are 0 <= x < 1
    aindex = np.where(((spos[:, 0] >= 0.0) & (spos[:, 0] < 1.0)
                     & (spos[:, 1] >= 0.0) & (spos[:, 1] < 1.0)
                     & (spos[:, 2] >= 0.0) & (spos[:, 2] < 1.0)))
    assert len(aindex[0]) == newnatoms, 'Filtering failed: ' + str(newnatoms) + 'atoms expected, ' + str(len(aindex[0])) + ' found'
    
    # Build new system
    newsystem = atomman.core.System(atoms=system2.atoms[aindex], box=system2.box)
    
    # Normalize box vectors
    newsystem.box_set(a=newsystem.box.a, b=newsystem.box.b, c=newsystem.box.c,
                      alpha=newsystem.box.alpha, beta=newsystem.box.beta,
                      gamma=newsystem.box.gamma,
                      scale=True)
        
    # Wrap to ensure all atoms are within the box
    newsystem.wrap()
    
    return newsystem