# coding: utf-8
# Standard Python libraries
from copy import deepcopy
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
import atomman.unitconvert as uc
from ... import Box, System
from ...region import PlaneSet

def array_boundary(self,
                   box: Box,
                   width: float) -> PlaneSet:
    """
    Constructs a shape associated with the boundary regions used by the
    periodicarray() generation method.  The returned shape will encompass
    all atoms except those within width distance of the non-periodic
    surface.

    Parameters
    ----------
    box : atomman.Box
        The box associated with the full (base) system.
    width : float
        The width of the boundary region

    Returns
    -------
    atomman.region.PlaneSet
        The Shape object excluding the boundary region
    """
    # Get plane shapes for the two non-periodic directions
    planes = [box.planes[self.cutindex], box.planes[self.cutindex + 3]]

    # Shift plane points by width in the plane normal directions
    for plane in planes:
        plane.point -= width * plane.normal

    # Create and return shape
    return PlaneSet(planes)

def periodicarray(self,
                  sizemults: Optional[tuple] = None,
                  amin: float = 0.0,
                  bmin: float = 0.0,
                  cmin: float = 0.0,
                  shift: Optional[npt.ArrayLike] = None,
                  shiftindex: Optional[int] = None,
                  shiftscale: bool = False,
                  center: Optional[npt.ArrayLike] = None,
                  centerscale: bool = False,
                  boundarywidth: float = 0.0,
                  boundaryscale: bool = False,
                  linear: bool = False,
                  cutoff: Optional[float] = None,
                  return_base_system: bool = False
                  ) -> Union[System, Tuple[System, System]]:
    """
    Constructs a dislocation monopole atomic configuration containing a
    single perfectly straight dislocation. The resulting system will be
    periodic along the box vector direction that corresponds to the
    dislocation's line direction, and non-periodic in the other two box
    vector directions.  Boundary atoms near the two free surfaces will be
    identified by changing their atype values making it easy to identify
    them later for assigning different boundary conditions.

    Parameters
    ----------
    sizemults : tuple, optional
        The size multipliers to use when generating the system.  Values
        are limited to being positive integers.  The multipliers for the
        two non-periodic directions must be even.  If not given, the
        default multipliers will be 2 for the non-periodic directions and
        1 for the periodic direction.
    amin : float, optional
        A minimum thickness to use for the a box vector direction of the
        final system.  Default value is 0.0.  For the non-periodic
        directions, the resulting vector multiplier will be even.  If both
        amin and sizemults is given, then the larger multiplier for the two
        will be used.
    bmin : float, optional
        A minimum thickness to use for the b box vector direction of the
        final system.  Default value is 0.0.  For the non-periodic
        directions, the resulting vector multiplier will be even.  If both
        bmin and sizemults is given, then the larger multiplier for the two
        will be used.
    cmin : float, optional
        A minimum thickness to use for the c box vector direction of the
        final system.  Default value is 0.0.  For the non-periodic
        directions, the resulting vector multiplier will be even.  If both
        cmin and sizemults is given, then the larger multiplier for the two
        will be used.
    shift : float, optional
        A rigid body shift to apply to the rotated cell prior to inserting
        the dislocation.  Should be selected such that the ideal slip plane
        does not correspond to any atomic planes.  Is taken as absolute if
        shiftscale is False, or relative to the rotated cell's box vectors
        if shiftscale is True.  Cannot be given with shiftindex.  If
        neither shift nor shiftindex is given then shiftindex = 0 is used.
    shiftindex : float, optional
        The index of the identified optimum shifts based on the rotated
        cell to use.  Different values allow for the selection of different
        atomic planes neighboring the slip plane.  Note that shiftindex
        values only apply shifts normal to the slip plane; best shifts for
        non-planar dislocations (like bcc screw) may also need a shift in
        the slip plane.  Cannot be given with shiftindex.  If neither shift
        nor shiftindex is given then shiftindex = 0 is used.
    shiftscale : bool, optional
        If False (default), a given shift value will be taken as absolute
        Cartesian.  If True, a given shift will be taken relative to the
        rotated cell's box vectors.
    center : array-like object or None, optional
        Indicates where the dislocation is positioned in the system relative
        to the default position at (0, 0) along the box vectors associated
        with the dislocation solution's m- and n-axes.
    centerscale : bool, optional
        If False (default), a given center value will be taken as absolute
        Cartesian.  If True, a given center will be taken relative to the
        rotated cell's box vectors.
    boundarywidth : float, optional
        The width of the boundary region to apply.  Default value is 0.0,
        i.e. no boundary region.  All atoms in the boundary region will
        have their atype values changed and will be displaced by linear
        displacements.
    boundaryscale : bool, optional
        If False (Default), the boundarywidth will be taken as absolute.
        If True, the boundarywidth will be taken relative to the magnitude
        of the unit cell's a box vector.
    linear : bool, optional
        If True, then only linear displacements will be used and not the
        dislocation solution.  Using only linear displacements is useful
        for screw dislocations and dislocations with large stacking fault
        distances.  If False (default) then the dislocation solution will
        be used for the middle displacements and linear displacements only
        in the boundary region.
    cutoff : float, optional
        Cutoff distance to use for identifying duplicate atoms to remove.
        For dislocations with an edge component, applying the displacements
        creates an extra half-plane of atoms that will have (nearly)
        identical positions with other atoms after altering the boundary
        conditions.  Default value is 0.5 Angstrom.
    return_base_system : bool, optional
        If True then the dislocation-free base system corresponding to the
        dislocation system will also be returned.  The base system is used
        as a reference state for most of the dislocation analysis tools.

    Returns
    -------
    base_system : atomman.System
        The base "perfect crystal" reference system associated with the
        dislocation system.  If the Burgers vector has an edge component
        then the atoms deleted when generating disl_system will also be 
        deleted from base_system. Only returned if return_base_system is
        True.
    disl_system : atomman.System
        The generated periodic array of dislocations system.
    """
    # Set default sizemults
    if sizemults is None:
        sizemults = [2,2,2]
        sizemults[self.lineindex] = 1
    else:
        try:
            assert len(sizemults) == 3
            assert isinstance(sizemults[0], int) and sizemults[0] > 0
            assert isinstance(sizemults[1], int) and sizemults[1] > 0
            assert isinstance(sizemults[2], int) and sizemults[2] > 0
            assert sizemults[self.lineindex - 1] % 2 == 0
            assert sizemults[self.lineindex - 2] % 2 == 0
        except AssertionError as err:
            raise TypeError('Invalid sizemults: must be 3 positive integers, and the two not along the dislocation line must be even') from err

    # Adjust multipliers based on min parameters
    if amin > 0.0:
        amult = int(np.ceil(amin / self.rcell.box.a))
        if self.lineindex != 0 and amult % 2 == 1:
            amult += 1
        if amult > sizemults[0]:
            sizemults[0] = amult

    if bmin > 0.0:
        bmult = int(np.ceil(bmin / self.rcell.box.b))
        if self.lineindex != 1 and bmult % 2 == 1:
            bmult += 1
        if bmult > sizemults[1]:
            sizemults[1] = bmult

    if cmin > 0.0:
        cmult = int(np.ceil(cmin / self.rcell.box.c))
        if self.lineindex != 2 and cmult % 2 == 1:
            cmult += 1
        if cmult > sizemults[2]:
            sizemults[2] = cmult

    # Modify the non-periodic size multipliers
    sizemults[self.lineindex] = (0, sizemults[self.lineindex])
    sizemults[self.lineindex - 1] = (-sizemults[self.lineindex - 1] // 2,
                                        sizemults[self.lineindex - 1] // 2)
    sizemults[self.lineindex - 2] = (-sizemults[self.lineindex - 2] // 2,
                                        sizemults[self.lineindex - 2] // 2)

    # Handle shift parameters
    if shift is not None or shiftindex is not None:
        self.set_shift(shift, shiftindex, shiftscale)
    shift = self.shift

    # Handle center parameter
    if center is None:
        center = np.array([0,0,0])
    else:
        center = np.asarray(center)
    if centerscale:
        center = self.rcell.box.vector_crystal_to_cartesian(center)

    # Handle boundary parameters
    if boundaryscale is True:
        boundarywidth = boundarywidth * self.ucell.box.a

    # Create the system where the dislocation will be inserted
    base_system = self.rcell.supersize(*sizemults)
    base_system.atoms.pos += shift
    base_system.wrap()

    # Create the dislocation system
    disl_system = self.build_disl_array(base_system, center, linear=linear,
                                        bwidth=boundarywidth, cutoff=cutoff)
    # Trim deleted atoms from base_system
    base_system = base_system.atoms_ix[disl_system.atoms.old_id]

    # Apply boundary region
    if boundarywidth > 0.0:

        shape = self.array_boundary(base_system.box, boundarywidth)

        # Change atypes of atoms outside box
        disl_system.atoms.atype[shape.outside(disl_system.atoms.pos)] += base_system.natypes
        disl_system.symbols = 2 * base_system.symbols

    self.set_systems(base_system, disl_system)

    if return_base_system:
        return base_system, disl_system
    else:
        return disl_system


def build_disl_array(self,
                     base_system: System,
                     center: npt.ArrayLike,
                     linear: bool = False,
                     bwidth: Optional[float] = None,
                     cutoff: Optional[float] = None,
                     ) -> System:
    """
    Method that converts a bulk crystal system into a periodic array of
    dislocations.  A single dislocation is inserted using a dislocation
    solution.  The system's box and pbc are altered such that the system is
    periodic and compatible across the two box vectors contained in the slip
    plane.  The third box vector is non-periodic, resulting in free surfaces
    parallel to the dislocation's slip plane.

    Parameters
    ----------
    base_system : atomman.System
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

    # Set default values
    if bwidth is None:
        bwidth = uc.set_in_units(10, 'angstrom')
    if cutoff is None:
        cutoff = uc.set_in_units(0.5, 'angstrom')

    # Extract dislocation solution values
    m = self.dislsol.m
    n = self.dislsol.n
    burgers = self.dislsol.burgers

    # Extract system values
    pos = base_system.atoms.pos
    vects = base_system.box.vects
    spos = base_system.atoms_prop(key='pos', scale=True)

    # Extract orientation values
    lineindex = self.lineindex
    cutindex = self.cutindex
    motionindex = self.motionindex

    # Check for atoms exactly on the slip plane
    if np.isclose(spos[:, cutindex], 0.5, rtol=0).sum() > 0:
        raise ValueError("atom positions found on slip plane: apply a coordinate shift")

    # ---------------------- Boundary modification -------------------------- #

    # Modify box vector in the motion direction by +- burgers/2
    newvects = deepcopy(vects)
    if burgers.dot(m) > 0:
        newvects[motionindex] -= burgers / 2
    else:
        newvects[motionindex] += burgers / 2
    newbox = Box(vects=newvects, origin=base_system.box.origin)

    # Make boundary condition perpendicular to slip plane non-periodic
    newpbc = [True, True, True]
    newpbc[cutindex] = False

    # Get length of system along motionindex
    # WHICH IS BEST!?!?!?
    length = np.abs(vects[motionindex].dot(m))
    #length = np.linalg.norm(newvects[motionindex])

    # -------------------- duplicate atom identification -------------------- #

    # Create test system to identify "duplicate" atoms
    testsystem = System(atoms=deepcopy(base_system.atoms), box=newbox,
                           pbc=newpbc, symbols=base_system.symbols)

    # Apply linear gradient shift to all atoms
    testsystem.atoms.pos += linear_displacement(pos - center, burgers, length, m, n)
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
    ii = np.ones(base_system.natoms, dtype=bool)
    ii[dup_atom_ids] = False

    # Count found duplicate atoms
    found = base_system.natoms - ii.sum()

    # Count expected number of duplicates based on volume change
    expected = base_system.natoms - (base_system.natoms * newbox.volume / base_system.box.volume)
    if np.isclose(expected, round(expected)):
        expected = int(round(expected))
    else:
        raise ValueError('expected number of atoms to delete not an integer: check burgers vector')

    # Compare found versus expected number of atoms
    if found != expected:
        raise ValueError('Deleted atom mismatch: expected %i, found %i. Adjust system dimensions and/or cutoff' %(expected, found))

    # ---------------------- Build dislocation system ----------------------- #

    # Generate new system with duplicate atoms removed
    disl_system = System(atoms=base_system.atoms[ii], box=newbox, pbc=newpbc,
                         symbols=base_system.symbols)

    # Define old_id so atoms in newsystem can be mapped back to system
    disl_system.atoms.old_id = np.where(ii)[0]

    if linear:
        # Use only linear displacements
        disp = linear_displacement(disl_system.atoms.pos - center, burgers, length, m, n)
    
    else:
        # Identify boundary atoms
        miny = base_system.box.origin.dot(n)
        maxy = miny + vects[cutindex].dot(n)
        if maxy < miny:
            miny, maxy = maxy, miny
        y = disl_system.atoms.pos.dot(n)
        ii = np.where((y <= miny + bwidth) | (y >= maxy - bwidth))
        
        # Use dislsol in middle and linear displacements at boundary
        disp = self.dislsol.displacement(disl_system.atoms.pos - center)
        disp[:, cutindex] -= disp[:, cutindex].mean()
        disp[ii] = linear_displacement(disl_system.atoms.pos[ii] - center, burgers,
                                       length, m, n)
    
    # Displace atoms and wrap
    disl_system.atoms.pos += disp
    disl_system.wrap()
    
    return disl_system
    
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
