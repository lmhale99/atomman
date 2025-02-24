# coding: utf-8
# Standard Python libraries
import io
from itertools import product
from typing import Optional, Tuple, Union

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from yabadaba.record import Record

from . import free_surface_basis
from ..tools import miller
from ..region import Plane
from .. import System, load
from ..library import load_record, Database

try:
    import spglib
    has_spglib = True
except ImportError:
    has_spglib = False

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
                 shift: Optional[npt.ArrayLike] = None,
                 shiftindex: Optional[int] = None,
                 shiftscale: bool = False,
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
            Indicates the space lattice setting of the given unit cell, i.e.
            'p' for primitive, 'i' for body-centered, 'f' for face-centered,
            'a', 'b', or 'c' for side-centered and 't1', or 't2' for trigonal
            in a hexagonal setting.  Setting this with the appropriate
            conventional unit cell allows for identifying lattice vectors that
            are not integers with respect to the conventional unit cell.  This
            also creates the rotated cell from a compatible primitive cell,
            thereby the final dislocation configurations can be smaller than
            possible solely from the conventional unit cell.
        shift : array-like object, optional
            Applies a shift to all atoms. Different values allow for free surfaces with
            different termination planes to be selected.  shift is taken as absolute
            if shiftscale is False, or relative to the rotated cell's box vectors
            if shiftscale is True.  Cannot be given with shiftindex.  If
            neither shift nor shiftindex is given then shiftindex = 0 is used.
        shiftindex : float, optional
            The index of the identified shifts based on the rotated
            cell to use.  Different values allow for the selection of different
            atomic planes neighboring the slip plane.  Cannot be given with shift.
            If neither shift nor shiftindex is given then shiftindex = 0 is used.
        shiftscale : bool, optional
            If False (default), a given shift value will be taken as absolute
            Cartesian.  If True, a given shift will be taken relative to the
            rotated cell's box vectors.
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
        _, unique_indices = np.unique(pos[:, cutindex].round(numdec), return_index=True)
        coords = pos[unique_indices, cutindex]

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
        self.__system = None
        self.__surfacearea = None

        # Set shift
        self.set_shift(shift=shift, shiftindex=shiftindex, shiftscale=shiftscale)

    @classmethod
    def fromrecord(cls,
                   record: Union[str, io.IOBase, DM, Record],
                   ucell: Union[System, str, io.IOBase],
                   maxindex: Optional[int] = None,
                   tol: float = 1e-7):
        """
        Construct a FreeSurface object based on parameters in a free_surface
        record and unit cell information.

        Parameters
        ----------
        record : atomman.library.record.FreeSurface, str, file-like object or DataModelDict
            A FreeSurface record object or the model contents for one.
        ucell : atomman.System
            The unit cell to use in generating the system.
        maxindex : int, optional
            Max uvw index value to use in identifying the best uvw set for the
            out-of-plane vector.  If not given, will use the largest absolute
            index between the given hkl and the initial in-plane vector guesses.
        tol : float, optional
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        # Create record object if needed
        if not isinstance(record, Record):
            record = load_record('free_surface', model=record)

        # Extract parameters in the record
        hkl = miller.fromstring(record.parameters['hkl'])
        shiftindex = int(record.parameters.get('shiftindex', 0))
        cutboxvector = record.parameters['cutboxvector']
        conventional_setting = record.parameters.get('conventional_setting', 'p')

        return cls(hkl=hkl, ucell=ucell, cutboxvector=cutboxvector,
                   maxindex=maxindex, shiftindex=shiftindex,
                   conventional_setting=conventional_setting, tol=tol)

    @classmethod
    def fromdatabase(cls,
                     name: Optional[str] = None,
                     ucell: Optional[System] = None,
                     database: Optional[Database] = None,
                     prompt: bool = True,
                     maxindex: Optional[int] = None,
                     tol: float = 1e-7,
                     **kwargs):
        """
        Construct a FreeSurface object based on record(s) retrieved from the
        reference database.

        Parameters
        ----------
        name : str or None, optional
            The name of the free_surface record to retrieve from the database.
            Alternatively, you can use any other query keyword arguments supported
            by the free_surface record style (see **kwargs below for more info).
        ucell : atomman.System or None, optional
            The unit cell to use in generating the system.  If None (default), then
            the crystal_prototype record that matches the defect's family setting
            will be loaded from the database.  Note that if None then the
            crystal-specific info (lattice constants and symbols) should be given
            here as kwargs (see below).
        database : atomman.library.Database or None, optional
            A Database object to use to fetch the records.  If None (default), then
            a new Database instance will be created.
        prompt : bool
            If prompt=True (default) then a screen input will ask for a selection
            if multiple matching free_surface (or crystal_prototype) records are
            found.  If prompt=False, then an error will be thrown if multiple
            matches are found.
        maxindex : int, optional
            Max uvw index value to use in identifying the best uvw set for the
            out-of-plane vector.  If not given, will use the largest absolute
            index between the given hkl and the initial in-plane vector guesses.
        tol : float, optional
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        **kwargs : any
            The recognized kwargs include the query keywords for free_surface
            records (key, id, family, hkl, shiftindex, cutboxvector), and the
            crystal-specific parameters recognized by the prototype load style
            (a, b, c, alpha, beta, gamma, symbols).  The non-trivial
            crystal-specific parameters should be for the crystal if ucell is
            not given above as the crystal prototype lacks this information.
        """
        # Initialize a database if needed
        if database is None:
            database = Database()

        # Extract ucell parameters from kwargs
        prototype_kwargs = {}
        prototype_kwargs_names = ['a', 'b', 'c', 'alpha', 'beta', 'gamma', 'symbols']
        for prototype_kwargs_name in prototype_kwargs_names:
            if prototype_kwargs_name in kwargs:
                prototype_kwargs[prototype_kwargs_name] = kwargs.pop(prototype_kwargs_name)

        # Fetch matching defect record
        record = database.get_record('free_surface', name=name, prompt=prompt, **kwargs)

        # Fetch crystal prototype unit cell if needed
        if ucell is None:
            ucell = load('prototype', name=record.family, **prototype_kwargs)
        else:
            if len(prototype_kwargs) > 0:
                raise ValueError('crystal-specific kwargs cannot be given with ucell')

        return cls.fromrecord(record=record, ucell=ucell, maxindex=maxindex, tol=tol)

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
        """
        numpy.ndarray : The conventional Miller or Miller-Bravais crystal
        vectors associated with the rcell box vectors.
        """
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
    def shift(self) -> np.ndarray:
        """
        numpy.NDArray : The particular shift value that will be or was used to
        construct the defect system
        """
        return self.__shift

    @property
    def system(self) -> System:
        """atomman.System : The built free surface system."""
        if self.__system is not None:
            return self.__system
        else:
            raise AttributeError('system not yet built. Use build_system() or surface().')

    @property
    def surfacearea(self) -> float:
        """float : The surface area of one of the hkl planes."""
        if self.__surfacearea is not None:
            return self.__surfacearea
        else:
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
                shiftindex: Optional[int] = None,
                shiftscale: bool = None,
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
            different termination planes to be selected.  shift is taken as absolute
            if shiftscale is False, or relative to the rotated cell's box vectors
            if shiftscale is True.  Cannot be given with shiftindex.  If
            neither shift nor shiftindex is given then the current value set to the
            shift attribute will be used.
        shiftindex : float, optional
            The index of the identified shifts based on the rotated
            cell to use.  Different values allow for the selection of different
            atomic planes neighboring the slip plane.  Cannot be given with shift.
            If neither shift nor shiftindex is given then the current value set to
            the shift attribute will be used.
        shiftscale : bool, optional
            If False (default), a given shift value will be taken as absolute
            Cartesian.  If True, a given shift will be taken relative to the
            rotated cell's box vectors.
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
        if shift is not None or shiftindex is not None:
            self.set_shift(shift=shift, shiftindex=shiftindex, shiftscale=shiftscale)
        shift = self.shift

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
                      atol: float = 1e-8,
                      return_indices: bool = False
                      ) -> Union[np.ndarray, Tuple[np.ndarray, list]]:
        """
        Use crystal symmetry operations to filter the list of shift values to
        only those that are symmetrically unique.  Note that the identified
        unique shifts can still result in the creation of energetically
        equivalent free surfaces if the free surface introduces a symmetry
        operation not present in the bulk crystal.

        Parameters
        ----------
        symprec: float
            The symmetry precision tolerance value used in spglib.
        trial_image_range: int, more than or equal to 1.
            Maximum cell images searched in finding translationally equivalent
            planes.  The default value is one, which corresponds to search the
            27 neighbor images, [-1, 1]^3.  The default value may not be
            sufficient for largely distorted lattice.
        atol: float
            The absolute tolerance used in comparing two crystal planes.
        return_indices: bool
            If True then the indices of shift that correspond to the
            identified unique shifts will be returned as well.  Default value
            is False (only return the shift vectors).

        Returns
        -------
        unique_shifts: np.ndarray, (# of unique shifts, 3)
            The symmetrically unique shift vectors.
        unique_indices: list
            The indices of shifts that correspond to the identified unique
            shifts.
        """
        if not has_spglib:
            raise ImportError("FreeSurface.unique_shifts requires spglib. Use `pip install spglib`")
        if (not isinstance(trial_image_range, int)) or (trial_image_range <= 0):
            raise ValueError("trial_image_range should be positive integer.")

        # Planes for each shift
        normal = np.zeros(3)
        normal[self.cutindex] = 1.0
        planes = [Plane(normal, shift) for shift in self.shifts]

        # Get symmetry operations of rotated ucell
        vects, positions, numbers = self.ucell.dump('spglib_cell')
        rotated_vects = np.dot(vects, self.transform.T)
        dataset = spglib.get_symmetry_dataset((rotated_vects, positions, numbers), symprec=symprec)

        # Compatibility fix for newspglib
        if hasattr(dataset, 'rotations'):
            dataset = {'rotations': dataset.rotations,
                       'translations': dataset.translations,
                       'primitive_lattice': dataset.primitive_lattice}

        # Convert operations to Cartesian
        operations = []
        vects_tinv = np.linalg.inv(rotated_vects.T)
        for rotation, translation in zip(dataset['rotations'], dataset['translations']):
            rotation_cart = np.dot(np.dot(rotated_vects.T, rotation), vects_tinv)
            translation_cart = np.inner(translation, rotated_vects.T)

            # It is sufficient to consider only symmetry operation that preserve the normal vector.
            if np.allclose(np.dot(rotation_cart, normal), normal):
                operations.append((rotation_cart, translation_cart))

        unique_shifts = []
        unique_indices = []

        # Use primitive vectors for the search if available
        try:
            primitive_vects = dataset['primitive_lattice']
        except KeyError:
            primitive_vects = rotated_vects

        # List of trial displacements to search for a translation between two planes
        # The range [-1, 1] may not be sufficient for largely distorted lattice.
        trial_images = list(product(range(-trial_image_range, trial_image_range + 1), repeat=3))

        def is_equivalent_by_primitive_vects(plane1, plane2):
            # Check if two planes can be transformed to each other by lattice vectors
            for image in trial_images:
                rotation = np.eye(3)
                translation = np.inner(image, primitive_vects.T)
                new_plane1 = plane1.operate(rotation, translation)
                if new_plane1.isclose(plane2, atol=atol):
                    return True
            return False

        for i in range(len(self.shifts)-1, -1, -1):
            plane_i = planes[i]

            # Obtain symmetrically equivalent planes with the i-th plane
            equivalent_planes = []
            for rotation, translation in operations:
                new_plane = plane_i.operate(rotation, translation)

                # If the new plane is already found, skip it.
                skip = False
                for old_plane in equivalent_planes:
                    if new_plane.isclose(old_plane, atol=atol):
                        skip = True
                        continue
                if skip:
                    continue
                equivalent_planes.append(new_plane)

            # Compare with remained shifts
            is_unique = True
            for j in range(i-1, -1, -1):
                plane_j = planes[j]
                # Here, the two planes have the normal vector.
                for plane in equivalent_planes:
                    if is_equivalent_by_primitive_vects(plane, plane_j):
                        is_unique = False
                        break

                if not is_unique:
                    break

            if is_unique:
                unique_indices.append(i)

        unique_indices.sort()
        unique_shifts = self.shifts[unique_indices]

        if return_indices is True:
            return unique_shifts, unique_indices
        return unique_shifts

    def set_shift(self,
                  shift: Optional[npt.ArrayLike] = None,
                  shiftindex: Optional[int] = None,
                  shiftscale: bool = False):
        """
        Directly set the shift value based on shiftindex, or shift and shiftscale.
        NOTE that the shift value can alternatively be set during class initialization
        or when surface() is called.

        Parameters
        ----------
        shift : array-like object, optional
            Applies a shift to all atoms. Different values allow for free surfaces with
            different termination planes to be selected.  shift is taken as absolute
            if shiftscale is False, or relative to the rotated cell's box vectors
            if shiftscale is True.  Cannot be given with shiftindex.  If
            neither shift nor shiftindex is given then shiftindex = 0 is used.
        shiftindex : float, optional
            The index of the identified shifts based on the rotated
            cell to use.  Different values allow for the selection of different
            atomic planes neighboring the slip plane.  Cannot be given with shift.
            If neither shift nor shiftindex is given then shiftindex = 0 is used.
        shiftscale : bool, optional
            If False (default), a given shift value will be taken as absolute
            Cartesian.  If True, a given shift will be taken relative to the
            rotated cell's box vectors.
        """
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
