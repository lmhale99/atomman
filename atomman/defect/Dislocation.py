# coding: utf-8
# Standard Python libraries
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

from . import (dislocation_system_basis, solve_volterra_dislocation,
               dislocation_array)
from ..tools import miller, vect_angle
from ..region import PlaneSet, Cylinder

def boolean(value):
    """
    Allows conversion of strings to Booleans.
    
    Parameters
    ----------
    value : str or bool
        If str, then 'true' and 't' become True and 'false' and 'f' become
        false. If bool, simply return the value.
        
    Returns
    -------
    bool
        Equivalent bool of value.
        
    Raises
    ------
    ValueError
        If value is unrecognized.
    """
    
    # Pass Boolean values through without changing
    if value is True:
        return True
    elif value is False:
        return False
    
    # Convert strings
    elif value.lower() in ['true', 't']:
        return True
    elif value.lower() in ['false', 'f']:
        return False
    
    # Issue error for invalid string
    else:
        raise ValueError('Invalid Boolean string')

class Dislocation():
    
    def __init__(self, ucell, C, burgers, ξ_uvw, slip_hkl, m=[0,1,0],
                 n=[0,0,1], shift=None, shiftindex=None, shiftscale=None, tol=1e-8):
        """
        Class initializer.  Solves the dislocation solution and rotates the
        given unit cell to the proper orientation.
        
        Parameters
        ----------
        ucell : atomman.System
            The unit cell to use as the seed for generating the dislocation
            monopole system.
        C : atomman.ElasticConstants
            The elastic constants associated with the bulk crystal structure
            for ucell.
        burgers : array-like object
            The dislocation's Burgers vector given as a Miller or
            Miller-Bravais vector relative to ucell.
        ξ_uvw : array-like object
            The dislocation's line direction given as a Miller or
            Miller-Bravais vector relative to ucell.
        slip_hkl : array-like object
            The dislocation's slip plane given as a Miller or Miller-Bravais
            plane relative to ucell.
        m : array-like object, optional
            The m unit vector for the dislocation solution.  m, n, and ξ
            (dislocation line) should be right-hand orthogonal.  Default value
            is [0,1,0] (y-axis).
        n : array-like object, optional
            The n unit vector for the dislocation solution.  m, n, and ξ
            (dislocation line) should be right-hand orthogonal.  Default value
            is [0,0,1] (z-axis). n is normal to the dislocation slip plane.
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
        tol : float
            A cutoff tolerance used with obtaining the dislocation solution.
            Only needs to be changed if there are issues with obtaining a
            solution.
        """
        # Generate the dislocation solution
        self.__dislsol = solve_volterra_dislocation(C, burgers, ξ_uvw=ξ_uvw,
                                                    slip_hkl=slip_hkl, m=m, n=n,
                                                    box=ucell.box, tol=tol)
        self.__transform = self.dislsol.transform
        
        # Compute the uvws and transformation matrix corresponding to the
        # dislocation solution
        self.__uvws = dislocation_system_basis(ξ_uvw, slip_hkl, m=m, n=n,
                                        box=ucell.box, tol=tol)
        
        # Rotate ucell to correspond to the dislocation solution
        self.__rcell = ucell.rotate(self.uvws)
        self.__ucell = ucell
        
        # Identify the lattice vector aligned with ξ
        a_angle = vect_angle(self.rcell.box.avect, self.dislsol.ξ)
        b_angle = vect_angle(self.rcell.box.bvect, self.dislsol.ξ)
        c_angle = vect_angle(self.rcell.box.cvect, self.dislsol.ξ)
        
        if np.isclose(a_angle, 0.0, atol=tol, rtol=0.0) or np.isclose(a_angle, 180.0, atol=tol, rtol=0.0):
            self.__lineindex = 0
        elif np.isclose(b_angle, 0.0, atol=tol, rtol=0.0) or np.isclose(b_angle, 180.0, atol=tol, rtol=0.0):
            self.__lineindex = 1
        elif np.isclose(c_angle, 0.0, atol=tol, rtol=0.0) or np.isclose(c_angle, 180.0, atol=tol, rtol=0.0):
            self.__lineindex = 2
        else:
            raise ValueError('A box vector of the rotated system must correspond with m X n = ξ')
        
        # Identify the lattice vector that is not in the slip plane
        adot = np.dot(self.rcell.box.avect, self.dislsol.n)
        bdot = np.dot(self.rcell.box.bvect, self.dislsol.n)
        cdot = np.dot(self.rcell.box.cvect, self.dislsol.n)
        
        if np.isclose(adot, 0.0, atol=tol, rtol=0.0) and np.isclose(bdot, 0.0, atol=tol, rtol=0.0):
            self.__cutindex = 2
        elif np.isclose(adot, 0.0, atol=tol, rtol=0.0) and np.isclose(cdot, 0.0, atol=tol, rtol=0.0):
            self.__cutindex = 1
        elif np.isclose(bdot, 0.0, atol=tol, rtol=0.0) and np.isclose(cdot, 0.0, atol=tol, rtol=0.0):
            self.__cutindex = 0
        else:
            raise ValueError('Only one box vector of the rotated system must have a component along n')
        
        # Define out of plane unit vector 
        ovect = np.zeros(3)
        ovect[self.cutindex] = 1.0
        
        # Get out of plane width
        rcellwidth = self.rcell.box.vects[self.cutindex, self.cutindex]
        
        # Get the unique coordinates normal to the plane
        pos = self.rcell.atoms.pos
        numdec = - int(np.floor(np.log10(tol)))
        coords = np.unique(pos[:, self.cutindex].round(numdec))

        # Add periodic replica if missing
        if not np.isclose(coords[-1] - coords[0], rcellwidth, rtol=0.0, atol=tol):
            coords = np.append(coords, coords[0] + rcellwidth)
            
        # Compute the shifts
        relshifts = rcellwidth - (coords[1:] + coords[:-1]) / 2
        relshifts[relshifts > rcellwidth] -= rcellwidth
        relshifts[relshifts < 0.0] += rcellwidth
        self.__shifts = np.outer(np.sort(relshifts), ovect)
        
        # Handle shift parameters
        if shift is not None:
            if shiftindex is not None:
                raise ValueError('shift and shiftindex cannot both be given')
            if shiftscale is True:
                self.__shift = miller.vector_crystal_to_cartesian(shift, self.rcell.box)
            else:
                self.__shift = np.asarray(shift)
                assert self.__shift.shape == (3,) 
        elif shiftindex is not None:
            self.__shift = self.shifts[shiftindex]
        else:
            self.__shift = self.shifts[0]

        # Set base_system and disl_system to None
        self.__base_system = None
        self.__disl_system = None

    @classmethod
    def fromref(cls, ucell, C, model, tol=1e-8):
        """
        Initializes a Dislocation object based on pre-defined dislocation
        parameters from a reference record.

        Parameters
        ----------
        ucell : atomman.System
            The unit cell to use as the seed for generating the dislocation
            monopole system.
        C : atomman.ElasticConstants
            The elastic constants associated with the bulk crystal structure
            for ucell.
        model : str, file-like object or DataModelDict
            The reference record containing dislocation parameters to use.
        tol : float
            A cutoff tolerance used with obtaining the dislocation solution.
            Only needs to be changed if there are issues with obtaining a
            solution.
        """
        # Load the record
        model = DM(model).find('dislocation')

        # Extract the dislocation parameters
        cp = model['calculation-parameter']
        slip_hkl = np.array(cp['slip_hkl'].split(), dtype=int)
        ξ_uvw = np.array(cp['ξ_uvw'].split(), dtype=int)
        burgers = np.array(cp['burgers'].split(), dtype=float)
        m = np.array(cp.get('m', '0 1 0').split(), dtype=float)
        n = np.array(cp.get('n', '0 0 1').split(), dtype=float)

        shift = cp.get('shift', None)
        shiftindex = cp.get('shiftindex', None)
        shiftscale = boolean(cp.get('shiftscale', False))
        
        if shift is not None:
            shift = np.array(shift.split(), dtype='float')
        elif shiftindex is not None:
            shiftindex = int(shiftindex)
            
        return cls(ucell, C, burgers, ξ_uvw, slip_hkl, m=m, n=n, shift=shift,
                   shiftindex=shiftindex, shiftscale=shiftscale, tol=tol)
    
    @property
    def dislsol(self):
        """atomman.defect.VolterraDislocation : The elastic dislocation solution"""
        return self.__dislsol
    
    @property
    def uvws(self):
        """numpy.NDArray : The 3x3 array of uvw Miller vectors used to rotate ucell to rcell"""
        return self.__uvws
    
    @property
    def transform(self):
        """numpy.NDArray : The 3x3 Cartesian transformation matrix associated with rotating from ucell to rcell"""
        return self.__transform
    
    @property
    def ucell(self):
        """atomman.System : The crystal unit cell used as the basis for constructing the dislocation system"""
        return self.__ucell
    
    @property
    def rcell(self):
        """atomman.System : The cell associated with rotating ucell to coincide with the dislocation solution"""
        return self.__rcell
    
    @property
    def lineindex(self):
        """int : The index of the box vector that coincides with the dislocation line: 0=a, 1=b, 2=c"""
        return self.__lineindex
    
    @property
    def cutindex(self):
        """int : The index of the box vector that is not within the slip plane: 0=a, 1=b, 2=c"""
        return self.__cutindex
    
    @property
    def shifts(self):
        """list : All identified shifts that will place the slip plane halfway between atomic planes"""
        return self.__shifts
    
    @property
    def shift(self):
        """numpy.NDArray : The particular shift value that will be or was used to construct the dislocation system"""
        return self.__shift

    @property
    def base_system(self):
        """atomman.System : The "perfect crystal" reference system associated with the dislocation system"""
        if self.__base_system is not None:
            return self.__base_system
        else:
            raise ValueError('base_system not built yet: must call monopole() or periodicarray() first')

    @property
    def disl_system(self):
        """atomman.System : The generated dislocation system"""
        if self.__disl_system is not None:
            return self.__disl_system
        else:
            raise ValueError('disl_system not built yet: must call monopole() or periodicarray() first')

    def box_boundary(self, box, width):
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
    
    def array_boundary(self, box, width):
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

    def cylinder_boundary(self, box, width):
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
    
    def monopole(self, sizemults=None, amin=0.0, bmin=0.0, cmin=0.0,
                 shift=None, shiftindex=None, shiftscale=False,
                 boundaryshape='cylinder', boundarywidth=0.0,
                 boundaryscale=False, return_base_system=False):
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
        shift : float, optional
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
            try:
                assert len(sizemults) == 3
                assert isinstance(sizemults[0], int) and sizemults[0] > 0
                assert isinstance(sizemults[1], int) and sizemults[1] > 0
                assert isinstance(sizemults[2], int) and sizemults[2] > 0
                assert sizemults[self.lineindex - 1] % 2 == 0
                assert sizemults[self.lineindex - 2] % 2 == 0
            except:
                raise TypeError('Invalid sizemults: must be 3 positive integers, and the two not along the dislocation line must be even')
        
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
        if shift is not None:
            if shiftindex is not None:
                raise ValueError('shift and shiftindex cannot both be given')
            if shiftscale is True:
                self.__shift = miller.vector_crystal_to_cartesian(shift, self.rcell.box)
            else:
                self.__shift = np.asarray(shift)
                assert shift.shape == (3,) 
        elif shiftindex is not None:
            self.__shift = self.shifts[shiftindex]
        shift = self.shift
        
        # Handle boundary parameters
        if boundaryscale is True:
            boundarywidth = boundarywidth * self.ucell.box.a
        if boundaryshape not in ['cylinder', 'box']:
            raise ValueError('boundaryshape must be "cylinder" or "box"')
        
        # Create the system where the dislocation will be inserted
        base_system = self.rcell.supersize(*sizemults)
        base_system.atoms.pos += shift
        base_system.wrap()
        self.__base_system = base_system
        
        # Copy the system and displace atoms according to the dislocation solution
        disl_system = deepcopy(base_system)
        disl_system.atoms.pos += self.dislsol.displacement(disl_system.atoms.pos)
        disl_system.pbc = [False, False, False]
        disl_system.pbc[self.lineindex] = True
        disl_system.wrap()
        self.__disl_system = disl_system
        
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

    def periodicarray(self, sizemults=None, amin=0.0, bmin=0.0, cmin=0.0,
                      shift=None, shiftindex=None, shiftscale=False,
                      boundarywidth=0.0, boundaryscale=False, linear=False,
                      cutoff=None, return_base_system=False):
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
            except:
                raise TypeError('Invalid sizemults: must be 3 positive integers, and the two not along the dislocation line must be even')
        
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
        if shift is not None:
            if shiftindex is not None:
                raise ValueError('shift and shiftindex cannot both be given')
            if shiftscale is True:
                self.__shift = miller.vector_crystal_to_cartesian(shift, self.rcell.box)
            else:
                self.__shift = np.asarray(shift)
                assert shift.shape == (3,) 
        elif shiftindex is not None:
            self.__shift = self.shifts[shiftindex]
        shift = self.shift
        
        # Handle boundary parameters
        if boundaryscale is True:
            boundarywidth = boundarywidth * self.ucell.box.a
            
        # Create the system where the dislocation will be inserted
        base_system = self.rcell.supersize(*sizemults)
        base_system.atoms.pos += shift
        base_system.wrap()
        
        # Create the dislocation system
        if linear is True:
            disl_system = dislocation_array(base_system, m=self.dislsol.m, n=self.dislsol.n,
                                            burgers=self.dislsol.burgers, cutoff=cutoff)
        else:
            disl_system = dislocation_array(base_system, self.dislsol,
                                            bwidth=boundarywidth, cutoff=cutoff)
        self.__disl_system = disl_system
        base_system = base_system.atoms_ix[disl_system.atoms.old_id]
        self.__base_system = base_system

        # Apply boundary region
        if boundarywidth > 0.0:
            
            shape = self.array_boundary(base_system.box, boundarywidth)
            
            # Change atypes of atoms outside box
            disl_system.atoms.atype[shape.outside(disl_system.atoms.pos)] += base_system.natypes
            disl_system.symbols = 2 * base_system.symbols
        
        if return_base_system:
            return base_system, disl_system
        else:
            return disl_system
        