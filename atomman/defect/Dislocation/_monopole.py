# coding: utf-8
# Standard Python libraries
from copy import deepcopy
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from ... import Box, System
from ...region import PlaneSet, Cylinder

def box_boundary(self,
                 box: Box,
                 width: float) -> PlaneSet:
    """
    Constructs a shape associated with the box-style boundary region.  Used
    by the monopole() generation method.  The returned shape will encompass
    all atoms except those within width distance of the two non-periodic
    surfaces.

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
    planes = []
    for i in range(3):
        if i == self.lineindex:
            continue
        planes.append(box.planes[i])
        planes.append(box.planes[i+3])

    # Shift plane points by width in the plane normal directions
    for plane in planes:
        plane.point -= width * plane.normal

    # Create and return shape
    return PlaneSet(planes)

def cylinder_boundary(self,
                      box: Box,
                      width: float) -> Cylinder:
    """
    Constructs a shape associated with the cylinder-style boundary region.
    Used by the monopole() generation method.  The returned shape will
    encompass a cylinder of atoms centered around the dislocation line
    leaving a boundary region that will be at least width wide everywhere.

    Parameters
    ----------
    box : atomman.Box
        The box associated with the full (base) system.
    width : float
        The minimum width of the boundary region.

    Returns
    -------
    atomman.region.Cylinder
        The Shape object excluding the boundary region
    """
    # Reduce the problem to 2D: solution independent of ξ position
    mn = np.array([self.dislsol.m, self.dislsol.n])
    vect1 = mn.dot(box.vects[self.lineindex - 2])
    vect2 = mn.dot(box.vects[self.lineindex - 1])
    origin = mn.dot(box.origin)

    # Compute normal vectors to box vectors
    normal_vect1 = np.array([-vect1[1], vect1[0]])
    normal_vect2 = np.array([vect2[1], -vect2[0]])
    normal_vect1 = normal_vect1 / np.linalg.norm(normal_vect1)
    normal_vect2 = normal_vect2 / np.linalg.norm(normal_vect2)

    def line(p1, p2):
        """
        Defines a 2D line as used by the intersection function
        """
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0]*p2[1] - p2[0]*p1[1])
        return A, B, -C

    def intersection(L1, L2):
        """
        Identifies the (x,y) coordinates where two 2D lines intersect
        """
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return x, y
        else:
            return False

    # Define normal lines as originating at (0,0)
    normal_line_1 = line([0,0], normal_vect1)
    normal_line_2 = line([0,0], normal_vect2)

    # Define boundary lines based on 2D box corners
    bound_bot1 = line(origin, origin + vect1)
    bound_bot2 = line(origin, origin + vect2)
    bound_top1 = line(origin + vect2, origin + vect2 + vect1)
    bound_top2 = line(origin + vect1, origin + vect1 + vect2)

    # Identify intersection points between normal lines and the boundary lines
    intersections = np.array([intersection(normal_line_1, bound_bot1),
                              intersection(normal_line_2, bound_bot2),
                              intersection(normal_line_1, bound_top1),
                              intersection(normal_line_2, bound_top2)])

    # Find distance between (0,0) and the closest intercept
    smallest = np.min(np.linalg.norm(intersections, axis=1))

    # Radius = smallest distance minus the boundary width
    radius =  smallest - width

    # Axis is along the line direction and includes point (0,0,0)
    center1 = np.zeros(3)
    center2 = box.vects[self.lineindex]

    return Cylinder(center1, center2, radius, endcaps=False)

def monopole(self,
             sizemults: Optional[tuple] = None,
             amin: float = 0.0,
             bmin: float = 0.0,
             cmin: float = 0.0,
             shift: Optional[npt.ArrayLike] = None,
             shiftindex: Optional[int] = None,
             shiftscale: bool = False,
             center: Optional[npt.ArrayLike] = None,
             centerscale: bool = False,
             boundaryshape: str = 'cylinder',
             boundarywidth: float = 0.0,
             boundaryscale: bool = False,
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
        The size multipliers to use when generating the system.  Values are
        limited to being positive integers.  The multipliers for the two
        non-periodic directions must be even.  If not given, the default
        multipliers will be 2 for the non-periodic directions and 1 for the
        periodic direction.
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
    shift : array-like object, optional
        A rigid body shift to apply to the rotated cell prior to inserting
        the dislocation.  Should be selected such that the ideal slip plane
        does not correspond to any atomic planes.  Is taken as absolute if
        shiftscale is False, or relative to the rotated cell's box vectors
        if shiftscale is True.  Cannot be given with shiftindex.  If
        neither shift nor shiftindex is given will use the shift set during
        class initialization.
    shiftindex : float, optional
        The index of the identified optimum shifts based on the rotated
        cell to use.  Different values allow for the selection of different
        atomic planes neighboring the slip plane.  Note that shiftindex
        values only apply shifts normal to the slip plane; best shifts for
        non-planar dislocations (like bcc screw) may also need a shift in
        the slip plane.  Cannot be given with shiftindex.  If neither shift
        nor shiftindex is given then shiftindex = 0 is used then will use
        the shift set during class initialization.
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
    boundaryshape : str, optional
        Indicates the shape of the boundary region to use.  Options are
        'cylinder' (default) and 'box'.  For 'cylinder', the non-boundary
        region is defined by a cylinder with axis along the dislocation
        line and a radius that ensures the boundary is at least
        boundarywidth thick.  For 'box', the boundary region will be
        exactly boundarywidth thick all around.
    boundarywidth : float, optional
        The width of the boundary region to apply.  Default value is 0.0,
        i.e. no boundary region.  All atoms in the boundary region will
        have their atype values changed.
    boundaryscale : bool, optional
        If False (Default), the boundarywidth will be taken as absolute.
        If True, the boundarywidth will be taken relative to the magnitude
        of the unit cell's a box vector.
    return_base_system : bool, optional
        If True then the dislocation-free base system corresponding to the
        dislocation system will also be returned.  The base system is used
        as a reference state for most of the dislocation analysis tools.

    Returns
    -------
    base_system : atomman.System
        The base "perfect crystal" reference system associated with the
        dislocation system. Only returned if return_base_system is True.
    disl_system : atomman.System
        The generated dislocation monopole system.
    """
    # Set default sizemults
    if sizemults is None:
        sizemults = [2,2,2]
        sizemults[self.lineindex] = 1
    else:
        sizemults = deepcopy(sizemults)
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
    if boundaryshape not in ['cylinder', 'box']:
        raise ValueError('boundaryshape must be "cylinder" or "box"')

    # Create the system where the dislocation will be inserted
    base_system = self.rcell.supersize(*sizemults)
    base_system.atoms.pos += shift
    base_system.wrap()

    # Copy the system and displace atoms according to the dislocation solution
    disl_system = deepcopy(base_system)
    disl_system.atoms.pos += self.dislsol.displacement(disl_system.atoms.pos - center)
    disl_system.pbc = [False, False, False]
    disl_system.pbc[self.lineindex] = True
    disl_system.wrap()

    self.set_systems(base_system, disl_system)

    # Apply boundary region
    if boundarywidth > 0.0:

        if boundaryshape == 'box':
            shape = self.box_boundary(base_system.box, boundarywidth)

        elif boundaryshape == 'cylinder':
            shape = self.cylinder_boundary(base_system.box, boundarywidth)

        # Change atypes of atoms outside box
        disl_system.atoms.atype[shape.outside(disl_system.atoms.pos)] += base_system.natypes
        disl_system.symbols = 2 * base_system.symbols

    if return_base_system:
        return base_system, disl_system
    else:
        return disl_system
