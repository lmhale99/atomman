# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy
import warnings

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..core import Atoms, Box, System
import atomman.unitconvert as uc

__all__ = ['dislocation_array']

def dislocation_array(system, dislsol=None, m=None, n=None, burgers=None, bwidth=None, cutoff=None):
    """
    Method that converts a bulk crystal system into a periodic array of
    dislocations.  A single dislocation is inserted using a dislocation
    solution.  The system's box and pbc are altered such that the system is
    periodic and compatible across the two box vectors contained in the slip
    plane.  The third box vector is non-periodic, resulting in free surfaces
    parallel to the dislocation's slip plane.
    
    Parameters
    ----------
    system : atomman.System
        A perfect, bulk atomic system.
    dislsol : atomman.defect.Stroh or atomman.defect.IsotropicVolterra, optional
        A dislocation solution to use to displace atoms by.  If not given,
        all atoms will be given linear displacements associated with the
        long-range limits.
    m : array-like object, optional
        The dislocation solution m unit vector.  This vector is in the slip
        plane and perpendicular to the dislocation line direction.  Only needed
        if dislsol is not given.
    n : array-like object, optional
        The dislocation solution n unit vector.  This vector is normal to the 
        slip plane.  Only needed if dislsol is not given.
    burgers : array-like object, optional
        The Cartesian Burger's vector for the dislocation relative to the
        given system's Cartesian coordinates.  Only needed if dislsol is not
        given.
    bwidth : float, optional
        The width of the boundary region at the free surfaces.  Atoms within
        the boundaries will be displaced by linear displacements instead of
        by the dislocation solution.  Only given if dislsol is not None.
        Default value if dislsol is given is 10 Angstroms.
    cutoff : float, optional
        Cutoff distance to use for identifying duplicate atoms to remove.
        For dislocations with an edge component, applying the displacements
        creates an extra half-plane of atoms that will have (nearly) identical
        positions with other atoms after altering the boundary conditions.
        Default cutoff value is 0.5 Angstrom.
    """
    # Input parameter setup
    if dislsol is None:
        if m is None or n is None:
            raise ValueError('m and n are needed if no dislsol is given')
        try:
            assert np.isclose(np.linalg.norm(m), 1.0)
            assert np.isclose(np.linalg.norm(n), 1.0)
            assert np.isclose(m.dot(n), 0.0)
        except:
            raise ValueError('m and n must be perpendicular unit vectors')
        if bwidth is not None:
            raise ValueError('bwidth not allowed if dislsol is not given')
    else:
        m = dislsol.m
        n = dislsol.n
        burgers = dislsol.burgers
        if bwidth is None:
            bwidth = uc.set_in_units(10, 'angstrom')
    if cutoff is None:
        cutoff = uc.set_in_units(0.5, 'angstrom')
    
    # Extract system values
    pos = system.atoms.pos
    vects = system.box.vects
    spos = system.atoms_prop(key='pos', scale=True)
    
    # Identify system orientation indices
    motionindex = None
    lineindex = None
    pnormindex = None
    for i in range(3):
        if np.isclose(np.abs(vects[i].dot(np.cross(m, n))), np.linalg.norm(vects[i])):
            lineindex = i
        elif not np.isclose(vects[i].dot(n), 0.0):
            if pnormindex is not None:
                raise ValueError("Multiple box vectors have components normal to dislocation solution's slip plane")
            pnormindex = i
    if lineindex is None:
        raise ValueError("No box vectors found parallel to dislocation solution's line vector")
    if pnormindex is None:
        raise ValueError("No box vectors have components normal to dislocation solution's slip plane")
    motionindex = 3 - (lineindex + pnormindex)
    
    # Compute new box vects and pbc
    newvects = deepcopy(system.box.vects)
    if burgers[motionindex] > 0:
        newvects[motionindex] -= burgers / 2
    else:
        newvects[motionindex] += burgers / 2
    newbox = Box(vects=newvects, origin=system.box.origin)
    newpbc = [True, True, True]
    newpbc[pnormindex] = False
    
    # Generate a test system to identify duplicate atoms
    length = np.abs(vects[motionindex].dot(m))
    testsystem = System(atoms=deepcopy(system.atoms), box=newbox, pbc=newpbc, symbols=system.symbols)
    testsystem.atoms.pos += linear_displacement(pos, burgers, length, m, n)
    testsystem.atoms.old_id = range(testsystem.natoms)
    
    # Identify boundary atoms to check
    spos = testsystem.atoms_prop(key='pos', scale=True)
    sburgers = 2 * burgers[motionindex] / (length)
    boundaryatoms = testsystem.atoms[(spos[:, motionindex] < sburgers) | (spos[:, motionindex] > 1.0 - sburgers)]
    
    # Compare distances between boundary atoms to identify duplicates
    dup_atom_ids = []
    mins = []
    for ni, i in enumerate(boundaryatoms.old_id[:-1]):
        js = boundaryatoms.old_id[ni+1:]
        try:
            distances = np.linalg.norm(testsystem.dvect(i, js), axis=1)
            mindistance = distances.min()
        except:
            mindistance = np.linalg.norm(testsystem.dvect(i, js))
        if mindistance < cutoff:
            dup_atom_ids.append(i)
    ii = np.ones(system.natoms, dtype=bool)
    ii[dup_atom_ids] = False
    
    # Generate new system with duplicate atoms removed
    newsystem = System(atoms=system.atoms[ii], box=newbox, pbc=newpbc, symbols=system.symbols)
    
    # Check if number of atoms deleted matches expected value
    expected = len(pos[(pos[:, motionindex] >= 0.0) & (pos[:, motionindex] <= np.abs(burgers[motionindex]))]) // 2
    actual = system.natoms - newsystem.natoms
    if expected != actual:
        warnings.warn('%i deleted atoms expected but %i atoms deleted' %(expected, actual), Warning)
    
    if dislsol is None:
        # Use only linear displacements
        disp = linear_displacement(newsystem.atoms.pos, burgers, length, m, n)
    
    else:
        # Identify boundary atoms
        miny = system.box.origin.dot(n)
        maxy = miny + vects[pnormindex].dot(n)
        if maxy < miny:
            miny, maxy = maxy, miny
        y = newsystem.atoms.pos.dot(n)
        ii = np.where((y <= miny + bwidth) | (y >= maxy - bwidth))
        
        # Use dislsol in middle and linear displacements at boundary
        disp = dislsol.displacement(newsystem.atoms.pos)
        disp[:, pnormindex] -= disp[:, pnormindex].mean()
        disp[ii] = linear_displacement(newsystem.atoms.pos[ii], burgers, length, m, n)
    
    # Displace atoms and wrap
    newsystem.atoms.pos += disp
    newsystem.wrap()
    
    return newsystem
    
def linear_displacement(pos, burgers, length, m, n):
    """
    Computes linear displacements associated with a dislocation in a system.
    
    Parameters
    ----------
    pos : array-like object
        List of Cartesian atomic positions.
    burgers : array-like object
        The Cartesian Burgers vector
    length : float
        The total length of the system along the m direction.
    m : array-like object, optional
        The dislocation solution m unit vector.  This vector is in the slip
        plane and perpendicular to the dislocation line direction.  Only needed
        if dislsol is not given.
    n : array-like object, optional
        The dislocation solution n unit vector.  This vector is normal to the 
        slip plane.  Only needed if dislsol is not given.
    """
    return np.outer(0.5 - np.sign(pos.dot(n)) * ((pos.dot(m) / (2 * length)) + 0.25), burgers)