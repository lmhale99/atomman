# coding: utf-8
# Standard Python libraries
from itertools import product
from typing import Optional

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from . import free_surface_basis
from ..tools import miller
from ..region import Plane
from .. import System

try:
    from spglib import get_symmetry_dataset
    spglib_loaded = True
except ImportError:
    spglib_loaded = False


class FreeSurface():
    """
    Class for generating free surface atomic configurations using clean planar slices.
    """

    def __init__(self,
                 hkl: npt.ArrayLike,
                 ucell: System,
                 cutboxvector: str = 'c',
                 maxindex: Optional[int] = None,
                 conventional_setting: str = 'p',
                 tol: float = 1e-7):
        """
        Class initializer.  Identifies the proper rotations for the given hkl plane
        and cutboxvector, and creates the rotated cell.

        Parameters
        ----------
        hkl : array-like object
            The free surface plane to generate expressed in either 3 indices
            Miller (hkl) format or 4 indices Miller-Bravais (hkil) format.
        ucell : atomman.System
            The unit cell to use in generating the system.
        cutboxvector : str, optional
            Specifies which of the three box vectors corresponds to the
            out-of-plane vector.  Default value is c.
        maxindex : int, optional
            Max uvw index value to use in identifying the best uvw set for the
            out-of-plane vector.  If not given, will use the largest absolute
            index between the given hkl and the initial in-plane vector guesses.
        conventional_setting : str, optional
            Allows for rotations of a primitive unit cell to be determined from
            (hkl) indices specified relative to a conventional unit cell.  Allowed
            settings: 'p' for primitive (no conversion), 'f' for face-centered,
            'i' for body-centered, and 'a', 'b', or 'c' for side-centered.  Default
            behavior is to perform no conversion, i.e. take (hkl) relative to the
            given ucell.
        tol : float, optional
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """

        # Pass parameters to free_surface_basis to get the rotation uvws
        uvws = free_surface_basis(hkl, box=ucell.box, cutboxvector=cutboxvector, maxindex=maxindex,
                                  conventional_setting=conventional_setting)

        # Generate the rotated cell
        # rcell.box.vects == uvws @ ucell.box.vects @ transform.T
        rcell, transform = ucell.rotate(uvws, return_transform=True)

        # Transform uvws to the conventional cell representation
        if uvws.shape == (3,3):
            uvws = miller.vector_primitive_to_conventional(uvws, conventional_setting)

        # Set cutindex and rcellwidth based on cutboxvector
        if cutboxvector == 'a':
            if rcell.box.bvect[0] != 0.0 or rcell.box.cvect[0] != 0.0:
                raise ValueError("box bvect and cvect cannot have x component for cutboxvector='a'")
            cutindex = 0

        elif cutboxvector == 'b':
            if rcell.box.avect[1] != 0.0 or rcell.box.cvect[1] != 0.0:
                raise ValueError("box avect and cvect cannot have y component for cutboxvector='b'")
            cutindex = 1

        elif cutboxvector == 'c':
            if rcell.box.avect[2] != 0.0 or rcell.box.bvect[2] != 0.0:
                raise ValueError("box avect and bvect cannot have z component for cutboxvector='c'")
            cutindex = 2

        # Define out of plane unit vector
        ovect = np.zeros(3)
        ovect[cutindex] = 1.0

        # Get out of plane width
        rcellwidth = rcell.box.vects[cutindex, cutindex]

        # Get the unique coordinates normal to the plane
        pos = rcell.atoms.pos
        numdec = - int(np.floor(np.log10(tol)))
        coords = np.unique(pos[:, cutindex].round(numdec))

        # Add periodic replica if missing
        if not np.isclose(coords[-1] - coords[0], rcellwidth, rtol=0.0, atol=tol):
            coords = np.append(coords, coords[0] + rcellwidth)

        # Compute the shifts
        relshifts = rcellwidth - (coords[1:] + coords[:-1]) / 2
        relshifts[relshifts > rcellwidth] -= rcellwidth
        relshifts[relshifts < 0.0] += rcellwidth
        shifts = np.outer(np.sort(relshifts), ovect)

        # Save attributes
        self.__hkl = np.asarray(hkl)
        self.__ucell = ucell
        self.__rcell = rcell
        self.__cutboxvector = cutboxvector
        self.__cutindex = cutindex
        self.__uvws = uvws
        self.__rcellwidth = rcellwidth
        self.__shifts = shifts
        self.__transform = transform
        self.__conventional_setting = conventional_setting

    @property
    def hkl(self) -> np.ndarray:
        """numpy.ndarray : Crystal plane in Miller or Miller-Bravais indices"""
        return self.__hkl

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell to use in building the defect system."""
        return self.__ucell

    @property
    def rcell(self) -> System:
        """atomman.System : the rotated cell to use in building the defect system."""
        return self.__rcell

    @property
    def cutboxvector(self) -> str:
        """str : The box vector for the cut direction."""
        return self.__cutboxvector

    @property
    def cutindex(self) -> int:
        """int : The Cartesian index for the cut direction."""
        return self.__cutindex

    @property
    def uvws(self) -> np.ndarray:
        """numpy.ndarray : The conventional Miller or Miller-Bravais crystal vectors associated with the rcell box vectors."""
        return self.__uvws

    @property
    def rcellwidth(self) -> float:
        """float : The width of rcell in the cutindex direction."""
        return self.__rcellwidth

    @property
    def shifts(self) -> list:
        """list : All shift values that place the fault halfway between atomic layers in rcell."""
        return self.__shifts

    @property
    def system(self) -> System:
        """atomman.System : The built free surface system."""
        try:
            return self.__system
        except:
            raise AttributeError('system not yet built. Use build_system() or surface().')

    @property
    def surfacearea(self) -> float:
        """float : The surface area of one of the hkl planes."""
        try:
            return self.__surfacearea
        except:
            raise AttributeError('system not yet built. Use build_system() or surface().')

    @property
    def transform(self) -> np.ndarray:
        """numpy.ndarray : The Cartesian transformation tensor associated with rotating from ucell to rcell"""
        return self.__transform

    @property
    def conventional_setting(self) -> str:
        """str : The lattice setting/centering associated with the conventional cell (used if ucell is primitive)"""
        return self.__conventional_setting

    def surface(self,
                shift: Optional[npt.ArrayLike] = None,
                vacuumwidth: Optional[float] = None,
                minwidth: Optional[float] = None,
                sizemults: Optional[list] = None,
                even: bool = False) -> System:
        """
        Generates and returns a free surface atomic system.

        Parameters
        ----------
        shift : array-like object, optional
            Applies a shift to all atoms. Different values allow for free surfaces with
            different termination planes to be selected.
        vacuumwidth : float, optional
            If given, the free surface is created by modifying the system's box to insert
            a region of vacuum with this width. This is typically used for DFT calculations
            where it is computationally preferable to insert a vacuum region and keep all
            dimensions periodic.
        sizemults : list or tuple, optional
            The three System.supersize multipliers [a_mult, b_mult, c_mult] to use on the
            rotated cell to build the final system. Note that the cutboxvector sizemult
            must be an integer and not a tuple.  Default value is [1, 1, 1].
        minwidth : float, optional
            If given, the sizemult along the cutboxvector will be selected such that the
            width of the resulting final system in that direction will be at least this
            value. If both sizemults and minwidth are given, then the larger of the two
            in the cutboxvector direction will be used.
        even : bool, optional
            A True value means that the sizemult for cutboxvector will be made an even
            number by adding 1 if it is odd.  Default value is False.

        Returns
        -------
        atomman.System
            The free surface atomic system.
        """

        # Set default function values
        if shift is None:
            shift = np.zeros(3)
        else:
            shift = np.asarray(shift)

        if sizemults is None:
            sizemults = [1, 1, 1]

        # Handle minwidth
        if minwidth is not None:
            mult = int(np.ceil(minwidth / self.rcellwidth))

            sizemult = sizemults[self.cutindex]
            if mult > np.abs(sizemult):
                sizemults[self.cutindex] = np.sign(sizemult) * mult

        # Handle even
        if even and sizemults[self.cutindex] % 2 == 1:
            if sizemults[self.cutindex] > 0:
                sizemults[self.cutindex] += 1
            else:
                sizemults[self.cutindex] -= 1

        # Define out of plane unit vector
        ovect = np.zeros(3)
        ovect[self.cutindex] = 1.0

        # Supersize and shift the system
        system = self.rcell.supersize(*sizemults)
        system.atoms.pos += shift
        system.wrap()

        # Change system's pbc
        system.pbc = [True, True, True]
        system.pbc[self.cutindex] = False

        # Insert vacuumwidth
        if vacuumwidth is not None:
            if vacuumwidth < 0:
                raise ValueError('vacuumwidth must be positive')
            newvects = system.box.vects
            newvects[self.cutindex, self.cutindex] += vacuumwidth
            neworigin = system.box.origin - ovect * vacuumwidth / 2
            system.box_set(vects=newvects, origin=neworigin)

        # Compute surfacearea based on cutboxvector
        if self.cutboxvector == 'a':
            surfacearea = np.linalg.norm(np.cross(system.box.bvect, system.box.cvect))

        elif self.cutboxvector == 'b':
            surfacearea = np.linalg.norm(np.cross(system.box.avect, system.box.cvect))

        elif self.cutboxvector == 'c':
            surfacearea = np.linalg.norm(np.cross(system.box.avect, system.box.bvect))

        # Save attributes
        self.__system = system
        self.__surfacearea = surfacearea

        return self.system

    def unique_shifts(self,
                      symprec: float = 1e-5,
                      trial_image_range: int = 1,
                      rtol: float = 1e-5,
                      atol: float = 1e-8) -> np.ndarray:
        """
        Return symmetrically nonequivalent shifts

        Parameters
        ----------
        symprec: float
            the tolerance value used in spglib
        trial_image_range: int, more than or equal to 1.
            Maximum cell images searched in finding translationally equivalent planes.
            The default value is one, which corresponds to search the 27 neighbor images, [-1, 1]^3.
            The default value may not be sufficient for largely distorted lattice.
        rtol: float
            the relative tolerance used in comparing two crystal planes
        atol: float
            the absolute tolerance used in comparing two crystal planes

        Returns
        -------
        unique_shifts: np.ndarray, (# of unique shifts, 3)
        """
        if not spglib_loaded:
            raise ImportError("FreeSurface.unique_shifts requires spglib. Use `pip install spglib`")
        if (not isinstance(trial_image_range, int)) or (trial_image_range <= 0):
            raise ValueError("trial_image_range should be positive integer.")

        # Get symmetry operations of rotated ucell
        lattice, positions, numbers = self.ucell.dump('spglib_cell')
        rotated_lattice = np.dot(lattice, self.transform.T)
        dataset = get_symmetry_dataset((rotated_lattice, positions, numbers), symprec=symprec)
        operations = []
        vects_tinv = np.linalg.inv(rotated_lattice.T)
        for rotation, translation in zip(dataset['rotations'], dataset['translations']):
            rotation_cart = np.dot(np.dot(rotated_lattice.T, rotation), vects_tinv)
            translation_cart = np.inner(translation, rotated_lattice.T)
            operations.append((rotation_cart, translation_cart))

        # Planes for each shift
        normal = np.zeros(3)
        normal[self.cutindex] = 1.0
        planes = [Plane(normal, shift) for shift in self.shifts]

        unique_shifts = []
        primitive_vects = dataset['primitive_lattice']

        # List of trial displacements to search for a translation between two planes
        # The range [-1, 1] may not be sufficient for largely distorted lattice.
        trial_images = list(product(range(-trial_image_range, trial_image_range + 1), repeat=3))

        def is_equivalent_by_primitive_vects(plane1, plane2):
            # Check if two planes can be transformed to each other by lattice vectors
            for image in trial_images:
                rotation = np.eye(3)
                translation = np.inner(image, primitive_vects.T)
                new_plane1 = plane1.operate(rotation, translation)
                if new_plane1.isclose(plane2, rtol=rtol, atol=atol):
                    return True
            return False

        for i in range(len(self.shifts)):
            plane_i = planes[i]

            # Obtain symmetrically equvalent planes with the i-th plane
            equivalent_planes = []
            for rotation, translation in operations:
                new_plane = plane_i.operate(rotation, translation)

                # new plane should preserve the normal vector
                if not np.allclose(new_plane.normal, normal):
                    continue

                # If the new plane is already found, skip it.
                if any([new_plane.isclose(plane, rtol=rtol, atol=atol) for plane in equivalent_planes]):
                    continue

                equivalent_planes.append(new_plane)

            # Compare with remained shifts
            is_unique = True
            for j in range(i + 1, len(self.shifts)):
                plane_j = planes[j]
                # Here, the two planes have the normal vector.
                for plane in equivalent_planes:
                    if is_equivalent_by_primitive_vects(plane, plane_j):
                        is_unique = False
                        break

                if not is_unique:
                    break

            if is_unique:
                unique_shifts.append(plane_i.point)

        unique_shifts = np.array(unique_shifts)
        return unique_shifts