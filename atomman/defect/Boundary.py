# coding: utf-8
# Standard Python imports
from typing import Optional, Tuple
from math import ceil

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from scipy.spatial.transform import Rotation

# Local imports
from .. import System
from ..tools import vect_angle, iaslist, miller


class Boundary():
    """
    Class for generating systems to investigate grain/phase boundaries
    """
    def __init__(self,
                 ucell1: System,
                 ucell2: System,
                 uvws1: npt.ArrayLike,
                 uvws2: npt.ArrayLike,
                 conventional_setting1: str = 'p',
                 conventional_setting2: str = 'p',
                 cutboxvector: str = 'c',
                 zerostrain: bool = False,
                 maxmult: int = 10,
                 tol: float = 1e-8):
        """
        Class initializer.  This is generic to allow for either a grain or
        phase boundary to be specified.

        Parameters
        ----------
        ucell1 : atomman.System
            The reference unit cell to use for the first grain.
        ucell2 : atomman.System
            The reference unit cell to use for the second grain.
        uvws1 : array-like object
            The three Miller(-Bravais) crystal vectors of ucell1 to use for
            orienting the first grain such that each crystal vector will be
            aligned with one of the box vectors of the final configuration.
        uvws2 : array-like object
            The three Miller(-Bravais) crystal vectors of ucell2 to use for
            orienting the second grain such that each crystal vector will be
            aligned with one of the box vectors of the final configuration.
        conventional_setting1 : str, optional
            Specifies which conventional lattice setting that ucell1 is in.
            The default value of 'p' takes ucell to be primitive, in which case
            the uvws1 values must be integers.  This must be specified
            in order to access non-integer lattice vectors for uvws1. 
        conventional_setting2 : str, optional
            Specifies which conventional lattice setting that ucell1 is in.
            The default value of 'p' takes ucell to be primitive, in which case
            the uvws1 values must be integers.  This must be specified
            in order to access non-integer lattice vectors for uvws1. 
        cutboxvector : str, optional
            Indicates which of the three box vectors that the boundary will
            be placed along. Default value is 'c'.
        zerostrain : bool, optional
            Setting this as True will check to see if the orientations of the
            two grains are fully compatible without introducing strain.  This
            is important for grain boundary energy comparisons.  Default value
            is False (no check performed).
        maxmult : int, optional
            The max integer multiplier to use for the zero strain check if
            zerostrain is True. 
        """
        # Set default values
        self.__mults1 = np.zeros(3)
        self.__mults2 = np.zeros(3)

        # Save given settings
        self.__ucell1 = ucell1
        self.__ucell2 = ucell2
        self.__uvws1 = uvws1
        self.__uvws2 = uvws2
        self.__cutboxvector = cutboxvector

        # Generate compatible primitive unit cells and associated transformation matrixes
        ucell_prim1, transform_c1_to_p1 = ucell1.dump('conventional_to_primitive',
                                                     setting=conventional_setting1,
                                                     return_transform=True, atol=tol)
        ucell_prim2, transform_c2_to_p2 = ucell2.dump('conventional_to_primitive',
                                                      setting=conventional_setting2,
                                                      return_transform=True, atol=tol)

        # Handle uvws inputs
        uvws1, uvws_prim1 = interpret_rotation_uvws(uvws1, conventional_setting=conventional_setting1)
        uvws2, uvws_prim2 = interpret_rotation_uvws(uvws2, conventional_setting=conventional_setting2)

        # Save primitive cell settings
        self.__ucell_prim1 = ucell_prim1
        self.__ucell_prim2 = ucell_prim2
        self.__transform_c1_to_p1 = Rotation.from_matrix(transform_c1_to_p1)
        self.__transform_c2_to_p2 = Rotation.from_matrix(transform_c2_to_p2)
        self.__uvws_prim1 = uvws_prim1
        self.__uvws_prim2 = uvws_prim2

        # Create rotated cells
        rcell1, transform_p1_to_r1 = ucell_prim1.rotate(uvws_prim1, return_transform=True)
        rcell2, transform_p2_to_r2 = ucell_prim2.rotate(uvws_prim2, return_transform=True)
        clean_wrap(rcell1)
        clean_wrap(rcell1)
        self.__rcell1 = rcell1
        self.__rcell2 = rcell2
        self.__transform_p1_to_r1 = Rotation.from_matrix(transform_p1_to_r1)
        self.__transform_p2_to_r2 = Rotation.from_matrix(transform_p2_to_r2)

        # Check boxvectors
        if cutboxvector == 'a':
            if rcell1.box.bvect[0] != 0.0 or rcell1.box.cvect[0] != 0.0:
                raise ValueError("box bvect and cvect cannot have x component for cutboxvector='a'")
            if rcell2.box.bvect[0] != 0.0 or rcell2.box.cvect[0] != 0.0:
                raise ValueError("box bvect and cvect cannot have x component for cutboxvector='a'")
            self.__cutindex = 0

        elif cutboxvector == 'b':
            if rcell1.box.avect[1] != 0.0 or rcell1.box.cvect[1] != 0.0:
                raise ValueError("box avect and cvect cannot have y component for cutboxvector='b'")
            if rcell2.box.avect[1] != 0.0 or rcell2.box.cvect[1] != 0.0:
                raise ValueError("box avect and cvect cannot have y component for cutboxvector='b'")
            self.__cutindex = 1

        elif cutboxvector == 'c':
            if rcell1.box.avect[2] != 0.0 or rcell1.box.bvect[2] != 0.0:
                raise ValueError("box avect and bvect cannot have z component for cutboxvector='c'")
            if rcell2.box.avect[2] != 0.0 or rcell2.box.bvect[2] != 0.0:
                raise ValueError("box avect and bvect cannot have z component for cutboxvector='c'")
            self.__cutindex = 2

        for i in range(3):
            if i == self.cutindex:
                continue
            if not np.isclose(vect_angle(rcell1.box.vects[i], rcell2.box.vects[i]), 0.0):
                raise ValueError('in-plane vectors of the two grains must be parallel')

        if zerostrain:
            strain = self.identifymults(maxmult)[2]
            if not np.allclose(strain, np.zeros(3)):
                raise ValueError('no zero strain configuration found') 

    @property
    def uvws1(self) -> np.ndarray:
        """numpy.NDArray: The three crystal vectors used to orient the first grain expressed as Miller(-Bravais) vectors of ucell1"""
        return self.__uvws1

    @property
    def uvws2(self) -> np.ndarray:
        """numpy.NDArray: The three crystal vectors used to orient the second grain expressed as Miller(-Bravais) vectors of ucell2"""
        return self.__uvws2

    @property
    def uvws_prim1(self) -> np.ndarray:
        """numpy.NDArray: The three crystal vectors used to orient the first grain expressed as Miller vectors of ucell_prim1"""
        return self.__uvws_prim1

    @property
    def uvws_prim2(self) -> np.ndarray:
        """numpy.NDArray: The three crystal vectors used to orient the second grain expressed as Miller vectors of ucell_prim2"""
        return self.__uvws_prim2

    @property
    def ucell1(self) -> System:
        """atomman.System: The conventional reference unit cell for the first grain."""
        return self.__ucell1

    @property
    def ucell2(self) -> System:
        """atomman.System: The conventional reference unit cell for the second grain."""
        return self.__ucell2
    
    @property
    def ucell_prim1(self) -> System:
        """atomman.System: The primitive reference unit cell for the first grain."""
        return self.__ucell_prim1

    @property
    def ucell_prim2(self) -> System:
        """atomman.System: The primitive reference unit cell for the second grain."""
        return self.__ucell_prim2

    @property
    def rcell1(self) -> System:
        """atomman.System: The rotated cell for the first grain."""
        return self.__rcell1

    @property
    def rcell2(self) -> System:
        """atomman.System: The rotated cell for the second grain."""
        return self.__rcell2

    @property
    def transform_c1_to_p1(self) -> Rotation:
        """
        scipy.spatial.transform.Rotation: The Cartesian rotation associated with ucell1 to ucell_prim1
        """
        return self.__transform_c1_to_p1

    @property
    def transform_c2_to_p2(self) -> Rotation:
        """
        scipy.spatial.transform.Rotation: The Cartesian rotation associated with ucell2 to ucell_prim2
        """
        return self.__transform_c2_to_p2

    @property
    def transform_p1_to_r1(self) -> Rotation:
        """
        scipy.spatial.transform.Rotation: The Cartesian rotation associated with ucell_prim1 to rcell1
        """
        return self.__transform_p1_to_r1

    @property
    def transform_p2_to_r2(self) -> Rotation:
        """
        scipy.spatial.transform.Rotation: The Cartesian rotation associated with ucell_prim2 to rcell2
        """
        return self.__transform_p2_to_r2

    @property
    def cutboxvector(self) -> str:
        """str: The box vector that the grain boundary is positioned along"""
        return self.__cutboxvector

    @property
    def cutindex(self) -> int:
        """int: The integer index associated with the cutboxvector setting"""
        return self.__cutindex

    @property
    def mults1(self) -> np.ndarray:
        """numpy.NDArray: Three int size multipliers for the top grain."""
        return self.__mults1

    @mults1.setter
    def mults1(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype=int)
        assert value.shape == (3,)
        self.__mults1 = value

    @property
    def mults2(self) -> np.ndarray:
        """numpy.NDArray: Three int size multipliers for the bottom grain."""
        return self.__mults2

    @mults2.setter
    def mults2(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype=int)
        assert value.shape == (3,)
        self.__mults2 = value

    def identifymults(self,
                      maxmult: int = 10,
                      minwidth: float = 0.0,
                      setvalues: bool = False
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compares the in-plane vectors of rcell1 and rcell2 to determine optimum
        system multipliers to minimize strain.

        Parameters
        ----------
        maxmult : int, optional
            The maximum size multiplier to use in searching for low strain
            states for the in-plane box vectors.  Both grains will be searched
            up to this max, so only one grain's multiplier is guaranteed to be
            maxmult or less.  Default value is 10.
        minwidth: float, optional
            A minimum width that is used to determine the size multipliers for
            the out-of-plane box vectors of the two grains.  This width limit is
            independently applied to both grains (not the sum) and is taken in
            the direction perpendicular to the grain boundary.  Default value
            is 0.0, which sets the multipliers to +-1.
        setvalues: bool, optional
            If True, the identified mults1 and mults2 values will be saved to
            the corresponding class attributes.  Default value is False.

        Returns
        -------
        mults1 : numpy.ndarray
            The suggested size multipliers to use with the top grain.
        mults2 : numpy.ndarray
            The suggested size multipliers to use with the bottom grain.
        strains : numpy.ndarray
            The estimated strains relative to the top grain that will result
            from the two grains positioned together with the suggested
            mults1 and mults2. Three values are given, one for each direction,
            with the out-of-plane strain always being reported as zero.
        """
        strains = np.zeros(3)
        mults1 = np.ones(3, dtype=int)
        mults2 = np.ones(3, dtype=int)

        for i in range(3):
            if i == self.cutindex:
                continue

            length1 = np.linalg.norm(self.rcell1.box.vects[i])
            length2 = np.linalg.norm(self.rcell2.box.vects[i])

            # If lengths are equal, both mults are 1 and strains are zero
            if length1 == length2:
                continue

            # Search for minimum strain configurations
            mult1 = 0
            mult2 = 0
            strain = 999999
            for multa in range(1, maxmult+1):

                # Check multiples of length1
                length = length1 * multa

                # Compute strain and compare to current known minimum
                multb = round(length / length2)
                teststrain = (length - multb * length2) / length
                if abs(teststrain) < abs(strain) and multb > 0:
                    mult1 = multa
                    mult2 = multb
                    strain = teststrain

                    if strain == 0.0:
                        break

                # Check multiples of length2
                length = length2 * multa

                # Compute strain and compare to current known minimum
                multb = round(length / length1)
                teststrain = (multb * length1 - length) / (multb * length1)

                if abs(teststrain) < abs(strain) and multb > 0:
                    mult1 = multb
                    mult2 = multa
                    strain = teststrain

                    if strain == 0.0:
                        break

            strains[i] = strain
            mults1[i] = mult1
            mults2[i] = mult2

        # Set out-of-plane multipliers based on minwidth
        if minwidth == 0.0:
            mults1[self.cutindex] = 1
            mults2[self.cutindex] = -1
        else:
            width1 = self.rcell1.box.vects[self.cutindex, self.cutindex]
            mults1[self.cutindex] = ceil(minwidth / width1)
            width2 = self.rcell1.box.vects[self.cutindex, self.cutindex]
            mults2[self.cutindex] = - ceil(minwidth / width2)

        if setvalues:
            self.mults1 = mults1
            self.mults2 = mults2

        return mults1, mults2, strains

    def boundary(self,
                 mults1: Optional[npt.ArrayLike] = None,
                 mults2: Optional[npt.ArrayLike] = None,
                 maxmult: Optional[int] = None,
                 minwidth: Optional[float] = None,
                 freesurface: bool = False,
                 straintype: str = 'top',
                 shift1: float = 0.0,
                 shift2: float = 0.0,
                 deleter: float = 0.1,
                 deletefrom: str = 'top'):
        """
        Generates a phase/grain boundary configuration.

        Parameters
        ----------
        mults1 : array-like object, optional
            Three int size multipliers to use with the top grain.  If given,
            mults2 is required and maxmults and minwidth cannot be given.  The
            mults1 class attribute will be updated to this value by the method.
            If mults1 is not given, then either maxmults must be given or the
            mults1, mults2 class attributes must be set prior to calling.
        mults2 : array-like object, optional
            Three int size multipliers to use with the bottom grain.  If given,
            mults1 is required and maxmults and minwidth cannot be given.  The
            mults2 class attribute will be updated to this value by the method.
            If mults2 is not given, then either maxmults must be given or the
            mults1, mults2 class attributes must be set prior to calling.
        maxmult : int, optional
            The maximum size multiplier to use in searching for low strain
            states for the in-plane box vectors.  Both grains will be searched
            up to this max, so only one grain's multiplier is guaranteed to be
            maxmult or less.  If given, mults1 and mults2 cannot be given.  If
            not given, then mults1 and mults2 must either be given or
            previously set as class attributes.
        minwidth: float, optional
            A minimum width that is used to determine the size multipliers for
            the out-of-plane box vectors of the two grains.  This width limit is
            independently applied to both grains (not the sum) and is taken in
            the direction perpendicular to the grain boundary.  Only allowed if
            maxmult is given.  If not given with maxmult, then will default to
            0.0.
        freesurface: bool, optional
            Indicates if the system is non-periodic along the cutboxvector and
            therefore contains a free surface (value of True) or if all box
            vectors are periodic and the system contains two boundaries (value
            of False).  Default value is False.
        straintype: str, optional
            Indicates how the lattice mismatch strain is applied.  'top' and
            'bottom' will strain only the associated grains while 'both' will
            divide the strain between both grains.
        shift1: float, optional
            A rigid body shift to apply along one of the two in-plane box
            vectors.  This is taken relative to the rcell's a vector or b
            vector if cutboxvector='a', so it should range between 0 and 1.
            Default value is 0.0.
        shift2: float, optional
            A rigid body shift to apply along one of the two in-plane box
            vectors.  This is taken relative to the rcell's c vector or b
            vector if cutboxvector='c', so it should range between 0 and 1.
            Default value is 0.0.
        deleter: float, optional
            Any atoms in the separate grains closer than this distance will
            have one atom of the pair deleted to prevent overlapping atoms.
            Varying this value can possibly result in a lower energy
            configuration.  Default value is 0.1.
        deletefrom: str, optional
            Indicates which grain 'top' or 'bottom' the atoms found with
            deletewidth will be deleted from.  Typically doesn't matter
            for symmetric tilt grain boundaries but does matter for
            asymmetric grain boundaries or phase boundaries.
        
        Returns
        -------
        atomman.System
            The grain/phase boundary atomic configuration
        natoms1 : int
            The number of atoms in the top grain.  This can be used to determine
            which atoms are in which grain as the first natoms1 atoms in the
            returned system will be in the top grain, and the rest in the bottom
            grain.
        """
        if mults1 is not None or mults2 is not None:
            assert mults1 is not None and mults2 is not None, (
                'mults1 and mults2 must be given together')
            assert maxmult is None and minwidth is None, (
                'maxmult and minwidth cannot be given with mults1, mults2')
            self.mults1 = mults1
            self.mults2 = mults2

        elif maxmult is not None:
            if minwidth is None:
                minwidth = 0.0
            self.identifymults(maxmult=maxmult, minwidth=minwidth, setvalues=True)

        elif minwidth is not None:
            raise ValueError('maxmult must be given with minwidth')

        # Create systems for each grain
        assert not np.any(np.isclose(self.mults1, np.zeros(3))), 'mults cannot be zero'
        assert not np.any(np.isclose(self.mults2, np.zeros(3))), 'mults cannot be zero'
        assert self.mults1[self.cutindex] > 0, 'out-of-plane mult for top grain must be positive'
        assert self.mults2[self.cutindex] < 0, 'out-of-plane mult for bottom grain must be negative'
        system1 = self.rcell1.supersize(*self.mults1)
        system2 = self.rcell2.supersize(*self.mults2)

        # Build system
        self.applystrain(system1, system2, straintype)
        system = self.mergesystems(system1, system2, freesurface)
        self.applyshift(system, shift1, shift2, system1.natoms)
        newsystem, natoms1 = self.deleteoverlaps(system, deleter, deletefrom, system1.natoms)

        return newsystem, natoms1

    def applystrain(self, system1, system2, straintype):
        """
        Adjust the in-plane box vectors to apply the necessary strain for lattice
        compatibility.
        """
        if straintype == 'top':
            # Strain in-plane vects of top grain to match bottom grain vects
            vects = system2.box.vects
            vects[self.cutindex] = system1.box.vects[self.cutindex]
            system1.box_set(vects=vects, origin=system1.box.origin, scale=True)

        elif straintype == 'bottom':
            # Strain in-plane vects of bottom grain to match top grain vects
            vects = system1.box.vects
            vects[self.cutindex] = system2.box.vects[self.cutindex]
            system2.box_set(vects=vects, origin=system2.box.origin, scale=True)

        elif straintype == 'both':
            # Average vects (works as in-plane vects are parallel)
            vects = (system1.box.vects + system2.box.vects) / 2

            # Change vects while retaining out-of-plane vects and origins.
            vects[self.cutindex] = system1.box.vects[self.cutindex]
            system1.box_set(vects=vects, origin=system1.box.origin, scale=True)
            vects[self.cutindex] = system2.box.vects[self.cutindex]
            system2.box_set(vects=vects, origin=system2.box.origin, scale=True)

        else:
            raise ValueError('Unknown straintype value: allowed values are top, bottom, or both')

    def mergesystems(self, system1, system2, freesurface):
        """
        """
        # Combine atoms
        system = system1.atoms_extend(system2.atoms)

        # Increase box size
        vects = system.box.vects
        vects[self.cutindex] += system2.box.vects[self.cutindex]
        system.box_set(vects=vects, origin=system2.box.origin)

        if freesurface:
            system.pbc[self.cutindex] = False

        return system

    def applyshift(self, system, shift1, shift2, natoms1):
        """
        """
        if self.cutboxvector == 'a':
            shift = shift1 * system.box.bvect + shift2 * system.box.cvect
        elif self.cutboxvector == 'b':
            shift = shift1 * system.box.avect + shift2 * system.box.cvect
        elif self.cutboxvector == 'c':
            shift = shift1 * system.box.avect + shift2 * system.box.bvect

        system.atoms.pos[:natoms1] += shift
        system.wrap()

    def deleteoverlaps(self, system, deleter, deletefrom, natoms1):
        """
        """
        if deleter == 0.0:
            return system, natoms1

        # Build neighborlist
        cutoff = 1.2
        if deleter > cutoff:
            cutoff = deleter + 0.2
        nlist = system.neighborlist(cutoff=cutoff)

        dup = set()
        for i in range(natoms1):
            if nlist.coord[i] == 0:
                continue

            dmag = system.dmag(i, nlist[i])
            if nlist.coord[i] == 1:
                if dmag < deleter:
                    if deletefrom == 'top':
                        dup.add(i)
                    elif deletefrom == 'bottom':
                        dup.update(nlist[i])
            else:
                dups = nlist[i][dmag < deleter]

                if deletefrom == 'top' and np.any(dups >= natoms1):
                    dup.add(i)
                elif deletefrom == 'bottom':
                    dup.update(dups[dups >= natoms1])

        keepindex = [x for i, x in enumerate(range(system.natoms)) if i not in dup]
        newsystem = system.atoms_ix[keepindex]

        if deletefrom == 'top':
            natoms1 = natoms1 - len(dup)

        return newsystem, natoms1

    def iterboundaryshift(self,
                          mults1: Optional[npt.ArrayLike] = None,
                          mults2: Optional[npt.ArrayLike] = None,
                          maxmult: Optional[int] = None,
                          minwidth: Optional[float] = None,
                          freesurface: bool = False,
                          straintype = 'top',
                          shifts1 = 0.0,
                          shifts2 = 0.0,
                          deleters = 0.1,
                          deletefrom = 'top'):
        """
        Generates multiple proposed phase/grain boundary configurations for the
        same boundary orientation where the in-plane shifts are varied as well
        as which overlapping atoms are deleted.

        Parameters
        ----------
        mults1 : array-like object, optional
            Three int size multipliers to use with the top grain.  If given,
            mults2 is required and maxmults and minwidth cannot be given.  The
            mults1 class attribute will be updated to this value by the method.
            If mults1 is not given, then either maxmults must be given or the
            mults1, mults2 class attributes must be set prior to calling.
        mults2 : array-like object, optional
            Three int size multipliers to use with the bottom grain.  If given,
            mults1 is required and maxmults and minwidth cannot be given.  The
            mults2 class attribute will be updated to this value by the method.
            If mults2 is not given, then either maxmults must be given or the
            mults1, mults2 class attributes must be set prior to calling.
        maxmult : int, optional
            The maximum size multiplier to use in searching for low strain
            states for the in-plane box vectors.  Both grains will be searched
            up to this max, so only one grain's multiplier is guaranteed to be
            maxmult or less.  If given, mults1 and mults2 cannot be given.  If
            not given, then mults1 and mults2 must either be given or
            previously set as class attributes.
        minwidth: float, optional
            A minimum width that is used to determine the size multipliers for
            the out-of-plane box vectors of the two grains.  This width limit is
            independently applied to both grains (not the sum) and is taken in
            the direction perpendicular to the grain boundary.  Only allowed if
            maxmult is given.  If not given with maxmult, then will default to
            0.0.
        freesurface: bool, optional
            Indicates if the system is non-periodic along the cutboxvector and
            therefore contains a free surface (value of True) or if all box
            vectors are periodic and the system contains two boundaries (value
            of False).  Default value is False.
        straintype: str, optional
            Indicates how the lattice mismatch strain is applied.  'top' and
            'bottom' will strain only the associated grains while 'both' will
            divide the strain between both grains.
        shifts1: int, float, or list, optional
            Indicates which rigid body shifts to apply along one of the two
            in-plane box vectors.  This is taken relative to the rcell's a
            vector or b vector if cutboxvector='a', so shift values should
            range between 0 and 1.  Giving a single float value will only
            perform one shift while giving a list of floats will iterate over
            all values.  If an int is given, then that number of equally-
            spaced shifts will be explored between 0 <= shift1 < 1.  Default
            value is 0.0 (no shift, no iteration).
        shifts2: int, float, or list, optional
            Indicates which rigid body shifts to apply along one of the two
            in-plane box vectors.  This is taken relative to the rcell's c
            vector or b vector if cutboxvector='c', so shift values should
            range between 0 and 1.  Giving a single float value will only
            perform one shift while giving a list of floats will iterate over
            all values.  If an int is given, then that number of equally-
            spaced shifts will be explored between 0 <= shift1 < 1.  Default
            value is 0.0 (no shift, no iteration).
        deleters: float or list, optional
            One or more interatomic spacing cutoffs to explore for
            identifying overlapping atoms between the two grains.  If multiple
            values are given, systems are yielded only if they differ from
            the previous system with the same in-plane shift.  Default value is
            0.1 (one value, no iteration).
        deletefrom: str, optional
            Indicates which grain 'top' or 'bottom' the atoms found with
            deletewidth will be deleted from.  A value of 'both' will iterate
            over both top and bottom deletions.  Typically doesn't matter
            for symmetric tilt grain boundaries but does matter for
            asymmetric grain boundaries or phase boundaries.
        
        Yields
        -------
        atomman.System
            The grain/phase boundary atomic configuration
        natoms1 : int
            The number of atoms in the top grain.  This can be used to determine
            which atoms are in which grain as the first natoms1 atoms in the
            returned system will be in the top grain, and the rest in the bottom
            grain.
        """

        # Set up values to iterate over
        shifts1 = self.interpret_shifts(shifts1)
        shifts2 = self.interpret_shifts(shifts2)
        deleters = [float(i) for i in iaslist(deleters)]
        if deletefrom == 'both':
            deletefroms = ['top', 'bottom']
        else:
            deletefroms = [deletefrom]

        # Iterate over all combos
        for shift1 in shifts1:
            for shift2 in shifts2:
                for deletefrom in deletefroms:

                    # deletewidth must be inner loop to do the natoms check
                    natoms = -999999
                    for deleter in deleters:
                        system, natoms1 = self.boundary(mults1=mults1,
                                                        mults2=mults2,
                                                        maxmult=maxmult,
                                                        minwidth=minwidth,
                                                        freesurface=freesurface,
                                                        straintype=straintype,
                                                        shift1=shift1,
                                                        shift2=shift2,
                                                        deleter=deleter,
                                                        deletefrom=deletefrom)

                        # Only yield systems where a different number of atoms was deleted
                        if system.natoms != natoms:
                            natoms = system.natoms
                            yield system, natoms1

    def interpret_shifts(self, val):
        """
        If int, do equal shifts in range [0, 1).
        If float or == 0, only do that value.
        If list, do those values
        """
        # Test for int value
        try:
            assert int(val) == val
            assert val > 0
        except (AssertionError, ValueError, TypeError):
            pass
        else:
            return np.linspace(0, 1, num=int(val), endpoint=False)

        return [float(i) for i in iaslist(val)]


def interpret_rotation_uvws(*uvws,
                            conventional_setting='p'):
    """
    Interprets the three rotation uvws for a crystal.
    
    Parameters
    ----------
    uvws : list, array
        The three Miller uvw vectors or Miller-Bravais uvtw vectors that the
        unit cell should be rotated to align with.  The three vectors can
        either be given directly or within a list.  Each value in the list can
        itself be an array-like object or a str representation of a Miller
        vector.
    conventional_setting : str, optional
        Indicates the conventional lattice setting associated with the Miller
        vectors given.  Default value is 'p' for a primitive unit cell.
        
    Returns
    -------
    uvws_conv : numpy.ndarray
        3x3 array of Miller vector indices relative to the conventional unit cell.
    uvws_prim : numpy.ndarray
        3x3 array of Miller vector indices relative to the primitive unit cell.
    
    Raises
    ------
    ValueError
        If the generated uvws_prim indices are not integer values.    
    """
    # Check for valid len of uvws
    if len(uvws) == 1:
        uvws = uvws[0] 
    if len(uvws) != 3:
        raise ValueError('invalid uvws: 3 values expected')
    
    # Convert given uvws into uvws_conv array
    uvws_conv = []
    for uvw in uvws:
        if isinstance(uvw, str):
            uvw = miller.fromstring(uvw)
        else:
            uvw = np.asarray(uvw, dtype=float)
            
        if uvw.shape == (4,):
            uvw = miller.vector4to3(uvw)
        elif uvw.shape != (3,):
            raise ValueError('invalid uvws: each uvw must have 3 or 4 indices')
        
        uvws_conv.append(uvw)
    uvws_conv = np.array(uvws_conv)
    
    
    # Convert uvws1 to the primitive cell
    uvws_prim = miller.vector_conventional_to_primitive(uvws_conv, setting=conventional_setting)
    
    # Check that uvws_prim are ints
    uvws_prim_int = np.array(np.rint(uvws_prim), dtype=int)
    if np.allclose(uvws_prim, uvws_prim_int):
        uvws_prim = uvws_prim_int
    else:
        raise ValueError('primitive Miller vector indices are not all ints, therefore they are not lattice vectors!')
    
    return uvws_conv, uvws_prim

def clean_wrap(system, atol=1e-7):
    """
    Wrap atoms around periodic boundaries cleanly

    Parameters
    ----------
    system : atomman.System
        The atomic system to wrap the atoms for.
    atol : float
        The absolute tolerance to use for identifying atoms on a boundary.
        Default value is 1e-7.
    """
    # Safely wrap atoms into the rcells
    system.atoms.pos += atol
    system.wrap()
    system.atoms.pos -= atol
    system.atoms.pos[np.isclose(system.atoms.pos, 0.0, atol=atol, rtol=0)] = 0
