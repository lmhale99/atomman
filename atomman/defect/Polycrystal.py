# Standard Python imports
from typing import Optional
from math import ceil
from copy import deepcopy

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

import scipy

# https://freud.readthedocs.io/en/latest/index.html
try:
    import freud
except:
    has_freud = False
else:
    has_freud = True

# atomman imports
from .. import load, Atoms, Box, System
from ..tools import random_rotation
from ..region import Plane, PlaneSet

class Grain():
    """
    Class that manages how individual grains in the polycrystal are defined
    and generated.
    """

    def __init__(self,
                 polycrystal,
                 point: Optional[npt.ArrayLike] = None,
                 rotation: Optional[npt.ArrayLike] = None,
                 scale: bool = False):
        """
        Defines a new grain with center point and rotation.

        Parameters
        ----------
        polycrystal : Polycrystal
            The parent polycrystal that the Grain is part of.
        point : array-like or None, optional
            The voronoi seed/reference position for the grain.  If not given,
            a random position within the polycrystal's box will be selected.
        rotation : array-like or None, optional
            A (3x3) rotation matrix defining the orientation for the grain's
            crystal.  If not given, a random rotation will be selected
            using an algorithm that uniformly samples the spherical rotation
            space.
        scale : bool, optional
            If False (default), point is taken as a Cartesian vector
            in the same reference as the polycrystal's box.  If True, point
            is taken as a relative (fractional) vector with respect to the
            polycrystal's box and will be scaled accordingly.
        """
        # Set polycrystal
        self.__polycrystal = polycrystal

        # Set point
        if point is None:
            point = self.polycrystal.rng.random(3)
            scale = True
        point = np.array(point)
        if point.shape != (3,):
            raise ValueError('point must be a 3D position vector')
        if scale:
            self.__point = self.polycrystal.box.position_relative_to_cartesian(point)
        else:
            self.__point = point

        # Set rotation
        if rotation is None:
            rotation = random_rotation(rng=self.polycrystal.rng)
        rotation = np.array(rotation)
        try:
            assert rotation.shape == (3,3)
            assert np.allclose(np.linalg.norm(rotation, axis=1), 1.0)
            assert np.allclose(np.inner(rotation[0], rotation[1:]), 0.0)
            assert np.allclose(np.cross(rotation[0], rotation[1]), rotation[2])
        except AssertionError as e:
            raise ValueError('invalid rotation matrix') from e
        self.__rotation = rotation

        # Set defaults
        self.__volume = None
        self.__region = None
        self.__extra_data = {}


    def solve_region_freud(self,
                           index: int,
                           voronoi,
                           freud_box):
        """
        Finds the plane set region for the grain based on voronoi output
        obtained from the freud package.
        
        Parameters
        ----------
        index : int
            The index number for the grain,.
        voronoi : freud.locality.Voronoi
            The Voronoi solution obtained from freud.
        freud_box : freud.box.Box
            The box represented in the freud Cartesian frame.  Used to
            change positions into the atomman Cartesian frame.
        """
        if not has_freud:
            raise ImportError('freud package required by this method: "pip install freud-analysis"')
        
        self.__volume = None
        self.__region = None
        self.__extra_data = {}

        point = self.point

        # Extract data for the grain from the voronoi solution
        freud_vertices = voronoi.polytopes[index]
        neighbor_dvects = voronoi.nlist.vectors[voronoi.nlist[:, 0] == index]
        self.__volume = voronoi.volumes[index]

        # Convert vertices and neighbor points to the atomman Cartesian frame
        vertices_system = load('freud', freud_box, freud_vertices,
                               origin=self.polycrystal.box.origin)
        vertices = vertices_system.atoms.pos
        neighbor_points = neighbor_dvects + point

        # Compute distance magnitude between point and vertices
        dmag_to_point = np.linalg.norm(vertices - point, axis=1)

        # Loop over neighboring seed positions
        planes = []
        for neigbor_point in neighbor_points:
            
            # Edge vertices are equidistant between the grain point and a neighbor point
            dmag_to_neighbor = np.linalg.norm(vertices - neigbor_point, axis=1)
            edge_verticies = vertices[np.isclose(dmag_to_point, dmag_to_neighbor)]
            if len(edge_verticies) < 3:
                continue

            # Define plane using the plane normal and a point in the plane
            plane_normal = np.cross(edge_verticies[1]-edge_verticies[0], edge_verticies[2]-edge_verticies[0])
            plane_point = edge_verticies[0]
            plane = Plane(normal=plane_normal, point=plane_point)
            
            # Flip sign on the normal if cell point is "above" the plane
            if plane.above(point):
                plane = Plane(normal=-plane_normal, point=plane_point)

            planes.append(plane)

        # Build the region using the list of planes
        self.__region = PlaneSet(planes=planes)

        # Save extra data for verification checks
        self.extra_data['vertices'] = vertices
        self.extra_data['freud_vertices'] = freud_vertices
        self.extra_data['neighbor_points'] = neighbor_points

    def solve_region_scipy(self,
                           index: int,
                           voronoi: scipy.spatial.Voronoi):
        """
        Finds the plane set region for the grain based on voronoi output
        obtained from the scipy package.
        
        Parameters
        ----------
        index : int
            The index number for the grain's point.  Note this will differ
            from the grain's index as it is based on a supercell system.
        voronoi : scipy.spatial.Voronoi
            The Voronoi solution obtained from scipy.
        """
        
        self.__volume = None
        self.__region = None
        self.__extra_data = {}

        # Loop over ridges
        planes = []
        all_edge_vertices = []
        for p_indices, v_indices in zip(voronoi.ridge_points, voronoi.ridge_vertices):
            
            # Check if ridge is associated with the point
            if index in p_indices:
                if -1 in v_indices:
                    v_indices = [x for x in v_indices if x != -1]
                if len(v_indices) < 3:
                    continue

                # Get the vertices positions in the ridge
                edge_verticies = voronoi.vertices[v_indices]
                all_edge_vertices.append(edge_verticies)

                # Define plane using the plane normal and a point in the plane
                plane_normal = np.cross(edge_verticies[1]-edge_verticies[0], edge_verticies[2]-edge_verticies[0])
                plane_point = edge_verticies[0]
                plane = Plane(normal=plane_normal, point=plane_point)
                
                # Flip sign on the normal if cell point is "above" the plane
                if plane.above(self.point):
                    plane = Plane(normal=-plane_normal, point=plane_point)

                planes.append(plane)
        self.__region = PlaneSet(planes)

        # Save extra data for verification checks
        self.extra_data['edge_vertices'] = all_edge_vertices

    def build_atoms(self,
                    refsystem: System) -> Atoms:
        """
        Builds atomic positions within the grain.  Requires that the region
        be defined.

        Parameters
        ----------
        refsystem : atomman.System
            The atomic system to select atoms from.  Should be large enough
            that the grain region will fit inside for any rotation.
        """
        if self.region is None:
            raise ValueError('grain region not solved yet!')
        # Rotate atoms in the reference system and shift them to center around seed point
        pos_rotated = refsystem.atoms.pos.dot(self.rotation)
        center_of_mass = pos_rotated.mean(axis=0)
        pos_rotated_shifted = pos_rotated - center_of_mass + self.point

        # Identify atoms within the grain 
        inside_grain = self.region.inside(pos_rotated_shifted)

        atoms = deepcopy(refsystem.atoms[inside_grain])
        atoms.pos = pos_rotated_shifted[inside_grain]

        return atoms

    @property
    def polycrystal(self):
        """Polycrystal : The parent polycrystal that the Grain is part of."""
        return self.__polycrystal

    @property
    def point(self) -> np.ndarray:
        """numpy.ndarray : The voronoi seed/reference position for the grain."""
        return self.__point
    
    @property
    def rotation(self) -> np.ndarray:
        """numpy.ndarray : The rotation matrix for orienting the grain."""
        return self.__rotation

    @property
    def region(self) -> Optional[PlaneSet]:
        """atomman.region.PlaneSet or None : Region object used for atom selection"""
        return self.__region
    
    @property
    def volume(self) -> Optional[float]:
        """float or None : The volume of the region"""
        return self.__volume

    @property
    def extra_data(self) -> dict:
        """dict: extra data fields for method verification checks"""
        return self.__extra_data

class Polycrystal():

    def __init__(self,
                 box: Box,
                 ucell: System,
                 ngrains: Optional[int] = None,
                 points: Optional[npt.ArrayLike] = None,
                 rotations: Optional[npt.ArrayLike] = None,
                 method: Optional[str] = None,
                 rng: Optional[np.random.Generator] = None,
                 seed: Optional[int] = None,
                 scale: bool = False
                 ):
        """
        Initialize an object that oversees the generation of polycrystal
        atomic configurations.

        Parameters
        ----------
        box : atomman.Box
            The cell box of the final atomic configuration to be generated.
        ucell : atomman.System
            The reference crystal unit cell to use for filling in the grain
            regions.
        ngrains : int or None, optional
            The number of grains to create.  Optional if points and/or
            rotations are given.
        points : array-like or None, optional
            An (Nx3) array-like object where N = ngrains that specifies the
            seed points used for the Voronoi analysis.  If not given, random
            points will be selected.
        rotations: array-like or None, optional
            An (Nx3x3) array-like object where N = ngrains that specifies the
            rotation transformation matrices to apply to ucell's positions for
            a given grain.  If not given, random rotations will be selected
            using an algorithm that uniformly samples the spherical rotation
            space.
        method : str or None, optional
            Specifies the method for the Voronoi calculation and grain region
            identification.  'freud' uses the freud package, 'scipy' uses the
            scipy package, and 'no' will only define the points and not solve.
            Default value of None uses 'freud' if the package is installed,
            and uses 'scipy' otherwise.
        rng : numpy.random.Generator or None, optional
            An existing random number generator to use.  If None
            (default), a new generator will be created based on seed.
        seed : int or None, optional
            A random number seed to use when creating a generator if rng is not
            given.  If None (default), then then fresh, unpredictable entropy
            will be pulled from the OS.
        scale : bool, optional
            Indicates if points (if given) are to be taken as Cartesian (False,
            default) or relative/fractional (True) with respect to box.
        """
        
        # Set required values
        if not isinstance(box, Box):
            raise TypeError('box must be an atomman.Box object')
        self.__box = box
        if not isinstance(ucell, System):
            raise TypeError('ucell must be an atomman.System object')
        self.__ucell = ucell

        # Manage initial rng, seed values
        if rng is None:
            self.__rng = np.random.default_rng(seed)
        elif seed is None:
            if not isinstance(rng, np.random.Generator):
                raise TypeError('rng must be a numpy.random.Generator object')
            self.__rng = rng
        else:
            raise ValueError('rng and seed cannot both be given')

        # Default values associated with an "unset" state
        self.__grains = []

        # Define the grains
        self.create_grains(ngrains=ngrains, points=points, rotations=rotations,
                           method=method, scale=scale)

    @property
    def box(self) -> Box:
        """atomman.Box : The cell box to use for the final system."""
        return self.__box
    
    @property
    def ucell(self) -> System:
        """atomman.System : The reference crystal unit cell to use for filling in the grain regions."""
        return self.__ucell
    
    @property
    def rng(self) -> np.random.Generator:
        """numpy.random.Generator :  Random number generator object"""
        return self.__rng

    @property
    def grains(self) -> list:
        """list : The separate Grain objects making up the polycrystal"""
        return self.__grains

    def create_grains(self,
                      box: Optional[Box] = None,
                      ngrains: Optional[int] = None,
                      points: Optional[npt.ArrayLike] = None,
                      rotations: Optional[npt.ArrayLike] = None,
                      method: Optional[str] = None,
                      rng: Optional[np.random.Generator] = None,
                      seed: Optional[int] = None,
                      scale: bool = False):
        """
        Defines the grains and rotations for the Voronoi polycrystal.
        This is automatically called during init, but can be called again
        later to generate new grains.

        Parameters
        ----------
        box : atomman.Box or None
            The cell box of the final atomic configuration to be generated.
            If None, will use the box set during init.
        ngrains : int or None, optional
            The number of grains to create.  If not given, will infer from
            the length of points or rotations.  If ngrains, points, and
            rotations are all not given, will use the number of grains
            previously used.
        points : array-like or None, optional
            An (Nx3) array-like object where N = ngrains that specifies the
            seed points used for the Voronoi analysis.  If not given, random
            points will be selected.
        rotations: array-like or None, optional
            An (Nx3x3) array-like object where N = ngrains that specifies the
            rotation transformation matrices to apply to ucell's positions for
            a given grain.  If not given, random rotations will be selected
            using an algorithm that uniformly samples the spherical rotation
            space.
        method : str or None, optional
            Specifies the method for the Voronoi calculation and grain region
            identification.  'freud' uses the freud package, 'scipy' uses the
            scipy package, and 'no' will only define the points and not solve.
            Default value of None uses 'freud' if the package is installed,
            and uses 'scipy' otherwise.
        rng : numpy.random.Generator or None, optional
            An existing random number generator to use.  If None
            (default), a new generator will be created based on seed.
        seed : int or None, optional
            A random number seed to use when creating a generator if rng is not
            given.  If None (default), then then fresh, unpredictable entropy
            will be pulled from the OS.
        scale : bool, optional
            Indicates if points (if given) are to be taken as Cartesian (False,
            default) or relative/fractional (True) with respect to box.
        """
        # Get/set box information
        if box is None:
            box = self.box
        elif not isinstance(box, Box):
            raise TypeError('box must be an atomman.Box object')
        else:
            self.__box = box

        # Get/set rng information
        if rng is None:
            if seed is None:
                rng = self.rng
            else:
                self.__rng = rng = np.random.default_rng(seed)
        else:
            if seed is not None:
                raise ValueError('rng and seed cannot both be given')
            elif not isinstance(rng, np.random.Generator):
                raise TypeError('rng must be a numpy.random.Generator object')
            else:
                self.__rng = rng

        # Check method values
        if method is None:
            if has_freud:
                method = 'freud'
            else:
                method = 'scipy'
        if method not in ['freud', 'scipy', 'no']:
            raise ValueError(f'unknown method style "{method}"')

        # Check dimensions of points and rotations
        if points is not None:
            points = np.array(points)
            if points.ndim != 2 or points.shape[1] != 3:
                raise ValueError('points must be an (Nx3) array')
        if rotations is not None:
            rotations = np.array(rotations)
            if rotations.ndim != 3 or rotations.shape[1] != 3 or rotations.shape[2] != 3:
                raise ValueError('rotations must be an (Nx3x3) array')
        
        # Check ngrains value
        if ngrains is None:
            if points is not None:
                # Get ngrains from number of points
                ngrains = points.shape[0]
            elif rotations is not None:
                # Get ngrains from number of rotations
                ngrains = rotations.shape[0]
            elif len(self.grains) > 0:
                # Get ngrains from previous number of grains
                ngrains = len(self.grains)
            
            if ngrains is None:
                raise ValueError('at least one of ngrains, points, and rotations must be provided')
            elif ngrains < 2:
                raise ValueError('need to have at least 2 grains')

        # Build None arrays if needed
        if points is None:
            points = np.array([None for i in range(ngrains)])
        if rotations is None:
            rotations = np.array([None for i in range(ngrains)])
        
        # Check that points and rotations both have ngrains values
        if points.shape[0] != ngrains or rotations.shape[0] != ngrains:
            raise ValueError('mismatch between number of grains, points and/or rotations')

        # Build new grains
        self.__grains = []
        for point, rotation in zip(points, rotations):
            self.grains.append(Grain(self, point=point, rotation=rotation, scale=scale))

        # Call a solve method
        if method == 'freud':
            self.solve_grains_freud()
        elif method == 'scipy':
            self.solve_grains_scipy()

    def solve_grains_freud(self):

        if not has_freud:
            raise ImportError('freud package required by this method: "pip install freud-analysis"')
        
        # Build an atomman system where the "atoms" are grain points
        points = np.empty((len(self.grains), 3))
        for i, grain in enumerate(self.grains):
            points[i,:] = grain.point
        points_system = System(atoms=Atoms(pos=points), box=self.box)
        
        # Convert atomman representation into the freud box, points representation
        freud_box, freud_points = points_system.dump('freud')

        # Initialize a voronoi object and solve based on the box and points
        voronoi = freud.locality.Voronoi()
        voronoi.compute((freud_box, freud_points))
        
        # Pass voronoi info to each grain to build planeset regions
        for i, grain in enumerate(self.grains):
            grain.solve_region_freud(i, voronoi, freud_box)

    def solve_grains_scipy(self):

        # Build an atomman system where the "atoms" are grain points
        points = np.empty((len(self.grains), 3))
        for i, grain in enumerate(self.grains):
            points[i,:] = grain.point
        points_system = System(atoms=Atoms(pos=points), box=self.box)
        
        # Create supercell to manage periodic boundaries
        points_supercell = points_system.supersize((-1,2),(-1,2),(-1,2))

        # Initialize a voronoi object and solve based on the supercell points
        voronoi = scipy.spatial.Voronoi(points_supercell.atoms.pos)

        # Find the voronoi point indices associated with the original points
        supercell_index = np.where(self.box.inside(voronoi.points))[0]
        
        # Pass voronoi info to each grain to build planeset regions
        for i, grain in zip(supercell_index, self.grains):
            grain.solve_region_scipy(i, voronoi)
            
            
    def build_polycrystal(self,
                          ucell: Optional[System] = None,
                          sizemults: Optional[tuple] = None,
                          verbose: bool = False) -> System:
        """
        Builds a polycrystal atomic configuration by filling in the grains
        and box with atoms.

        Parameters
        ----------
        ucell : atomman.System or None, optional
            Specifies the crystal unit cell to use for generating the grains.
            If None (default), will use the ucell set during class
            initialization.
        sizemults : tuple or None, optional
            System.supersize() size multipliers to use on ucell.  If given,
            the resulting system should be large enough that the grain region
            fits inside the supercell for any rotation.  Default value of None
            will guess sizemults for you.
        verbose : bool
            If set to True and the freud method was used, will print out
            the number of atoms in each grain compared to the estimated expected
            number based on ucell's density and the grain's volume.

        """

        if ucell is None:
            ucell = self.ucell

        density = ucell.natoms / ucell.box.volume

        # Find mults based on max distance between corners of box
        if sizemults is None:
            corners = np.array([
                self.box.origin,
                self.box.origin + self.box.avect,
                self.box.origin + self.box.bvect,
                self.box.origin + self.box.cvect,
                self.box.origin + self.box.avect + self.box.bvect,
                self.box.origin + self.box.avect + self.box.cvect,
                self.box.origin + self.box.bvect + self.box.cvect,
                self.box.origin + self.box.avect + self.box.bvect + self.box.cvect])
            maxwidth = np.linalg.norm(corners[1:] - corners[0], axis=1).max()
            sizemults = (ceil(maxwidth / ucell.box.a),
                         ceil(maxwidth / ucell.box.b),
                         ceil(maxwidth / ucell.box.c))
        
        # Generate reference system from ucell
        refsystem = ucell.supersize(*sizemults)
        
        # Build atoms for each grain
        all_atoms = None
        for grain in self.grains:
            atoms = grain.build_atoms(refsystem)
            if all_atoms is not None:
                all_atoms = all_atoms.extend(atoms)
            else:
                all_atoms = atoms

            if verbose and grain.volume is not None:
                natoms = atoms.natoms
                natoms_estimate = round(density * grain.volume)
                print(f'{natoms} atoms generated, {natoms_estimate} estimated')

        system = System(box=self.box, atoms=all_atoms, symbols=ucell.symbols)
        return system

            



