# coding: utf-8

# Standard Python libraries
from copy import deepcopy
from typing import Optional

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from ..core import Box, System
from . import VolterraDislocation
import atomman.unitconvert as uc

def dislocation_array(system: System,
                      dislsol: Optional[VolterraDislocation] = None,
                      m: Optional[npt.ArrayLike] = None,
                      n: Optional[npt.ArrayLike] = None,
                      burgers: Optional[npt.ArrayLike] = None,
                      bwidth: float = None,
                      cutoff: float = None) -> System:
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
    dislsol : atomman.defect.VolterraDislocation, optional
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
    
    Returns
    -------
    atomman.System
        The resulting periodic array of dislocations system.  An additional
        atoms property 'old_id' will be added to map the atoms in the defect
        system back to the associated atoms in the original system.
    """
    
    # ------------------------ Parameter handling --------------------------- #
    
    # Input parameter setup for linear gradients only
    if dislsol is None:
        if m is None or n is None:
            raise ValueError('m and n are needed if no dislsol is given')
        m = np.asarray(m)
        n = np.asarray(n)
        burgers = np.asarray(burgers)
        try:
            assert np.isclose(np.linalg.norm(m), 1.0)
            assert np.isclose(np.linalg.norm(n), 1.0)
            assert np.isclose(m.dot(n), 0.0)
        except:
            raise ValueError('m and n must be perpendicular unit vectors')
        if bwidth is not None:
            raise ValueError('bwidth not allowed if dislsol is not given')

    # Input parameter setup for dislsol
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
    
    # -------------------- Orientation identification ----------------------- #
    
    # Find box vector alligned with u = m x n
    u = np.cross(m,n)
    line = np.isclose(np.abs(np.dot(u, vects)), np.linalg.norm(vects, axis=1))
    if np.sum(line) != 1:
        raise ValueError('box vector aligned with u = m x n not found')
    lineindex = np.where(line)[0][0]

    # Find out-of-plane box vector
    out = ~np.isclose(np.dot(n, vects), 0.0)
    if np.sum(out) > 1:
        raise ValueError('multiple box vectors have out-of-plane components')
    pnormindex = np.where(out)[0][0]

    # Set third box vector as edge direction
    motionindex = 3 - (lineindex + pnormindex)

    # Check for atoms exactly on the slip plane
    if np.isclose(spos[:, pnormindex], 0.5, rtol=0).sum() > 0:
        raise ValueError("atom positions found on slip plane: apply a coordinate shift")
    
    # ---------------------- Boundary modification -------------------------- #
    
    # Modify box vector in the motion direction by +- burgers/2
    newvects = deepcopy(vects)
    if burgers.dot(m) > 0:
        newvects[motionindex] -= burgers / 2
    else:
        newvects[motionindex] += burgers / 2
    newbox = Box(vects=newvects, origin=system.box.origin)
        
    # Make boundary condition perpendicular to slip plane non-periodic
    newpbc = [True, True, True]
    newpbc[pnormindex] = False

    # Get length of system along motionindex 
    # WHICH IS BEST!?!?!?
    length = np.abs(vects[motionindex].dot(m))
    #length = np.linalg.norm(newvects[motionindex])
    
    # -------------------- duplicate atom identification -------------------- #
    
    # Create test system to identify "duplicate" atoms
    testsystem = System(atoms=deepcopy(system.atoms), box=newbox,
                           pbc=newpbc, symbols=system.symbols)

    # Apply linear gradient shift to all atoms
    testsystem.atoms.pos += linear_displacement(pos, burgers, length, m, n)
    testsystem.atoms.old_id = range(testsystem.natoms)

    # Identify atoms at the motionindex boundary to include in the duplicate check
    spos = testsystem.atoms_prop(key='pos', scale=True)
    sburgers = np.abs(2 * burgers[motionindex] / (length))
    boundaryatoms = testsystem.atoms[  (spos[:, motionindex] < sburgers) 
                                     | (spos[:, motionindex] > 1.0 - sburgers) ]

    # Compare distances between boundary atoms to identify duplicates
    dup_atom_ids = []
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

    # Count found duplicate atoms
    found = system.natoms - ii.sum()

    # Count expected number of duplicates based on volume change
    expected = system.natoms - (system.natoms * newbox.volume / system.box.volume)
    if np.isclose(expected, round(expected)):
        expected = int(round(expected))
    else:
        raise ValueError('expected number of atoms to delete not an integer: check burgers vector')
        
    # Compare found versus expected number of atoms
    if found != expected:
        raise ValueError('Deleted atom mismatch: expected %i, found %i. Adjust system dimensions and/or cutoff' %(expected, found))
    
    # ---------------------- Build dislocation system ----------------------- #
    
    # Generate new system with duplicate atoms removed
    newsystem = System(atoms=system.atoms[ii], box=newbox, pbc=newpbc,
                       symbols=system.symbols)
    
    # Define old_id so atoms in newsystem can be mapped back to system
    newsystem.atoms.old_id = np.where(ii)[0]

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
        disp[ii] = linear_displacement(newsystem.atoms.pos[ii], burgers,
                                       length, m, n)
    
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

    return np.outer(np.sign(pos.dot(n)) * (0.25 - pos.dot(m) / (2 * length)), burgers)