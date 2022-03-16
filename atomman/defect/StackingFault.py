# coding: utf-8
# Standard Python libraries
from copy import deepcopy
from typing import Generator, Optional, Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from ..tools import miller
from . import FreeSurface
from .. import System

class StackingFault(FreeSurface):
    """
    Class for generating stacking fault atomic configurations. 
    """
    
    def __init__(self,
                 hkl: npt.ArrayLike,
                 ucell: System,
                 cutboxvector: str = 'c',
                 maxindex: Optional[int] = None,
                 a1vect_uvw: Optional[npt.ArrayLike] = None,
                 a2vect_uvw: Optional[npt.ArrayLike] = None,
                 conventional_setting: str = 'p',
                 tol: float = 1e-8):
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
        a1vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.  If
            not given, will be set to the shortest in-plane lattice vector.
        a2vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.  If
            not given, will be set to the shortest in-plane lattice vector not
            parallel to a1vect_uvw.
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

        super().__init__(hkl, ucell, cutboxvector=cutboxvector, maxindex=maxindex,
                 conventional_setting=conventional_setting, tol=tol)
        
        # Extract a1vect and a2vect values from uvws and rcell.box
        if a1vect_uvw is None and a2vect_uvw is None:

            # Set a1index and a2index based on cutboxvector
            if self.cutboxvector == 'a':
                a1index, a2index = 1, 2
            elif self.cutboxvector == 'b':
                a1index, a2index = 2, 0
            elif self.cutboxvector == 'c':
                a1index, a2index = 0, 1

            self.a1vect_uvw = self.uvws[a1index]
            self.a2vect_uvw = self.uvws[a2index]
        
        # Set given a1vect_uvw and a2vect_uvw values
        elif a1vect_uvw is not None and a2vect_uvw is not None:
            self.a1vect_uvw = a1vect_uvw
            self.a2vect_uvw = a2vect_uvw
        
        else:
            raise ValueError('a1vect_uvw and a2vect_uvw either both need to be given or not given')
    
    @property
    def a1vect_uvw(self) -> np.ndarray:
        """numpy.ndarray : One of the two conventional lattice shift vectors in Miller or Miller-Bravais indices."""
        return self.__a1vect_uvw
    
    @property
    def a2vect_uvw(self) -> np.ndarray:
        """numpy.ndarray : One of the two conventional lattice shift vectors in Miller or Miller-Bravais indices."""
        return self.__a2vect_uvw
    
    @a1vect_uvw.setter
    def a1vect_uvw(self, value: npt.ArrayLike):
        
        value = np.asarray(value)

        # Check shape and convert [uvtw] to [uvw] if needed
        if value.shape == (4,):
            uvw = miller.vector4to3(value)
        elif value.shape == (3,):
            uvw = value
            if self.hkl.shape == (4,):
                value = miller.vector3to4(value)
        else:
            raise ValueError('Invalid uvw shape: must have 3 or 4 values.')

        # Convert to primitive cell
        uvw = miller.vector_conventional_to_primitive(uvw, self.conventional_setting)

        # Convert to Cartesian wrt rotated system
        cart = self.transform.dot(miller.vector_crystal_to_cartesian(uvw, self.ucell.box))

        # Check that Cartesian vector is in slip plane
        if not np.isclose(cart[self.cutindex], 0.0):
            raise ValueError(f'shift vector {value} not in fault plane {self.hkl}')

        # Save uvw and cart
        self.__a1vect_uvw = value
        self.__a1vect_cart = cart

    @a2vect_uvw.setter
    def a2vect_uvw(self, value: npt.ArrayLike):
        
        value = np.asarray(value)

        # Check shape and convert [uvtw] to [uvw] if needed
        if value.shape == (4,):
            uvw = miller.vector4to3(value)
        elif value.shape == (3,):
            uvw = value
            if self.hkl.shape == (4,):
                value = miller.vector3to4(value)
        else:
            raise ValueError('Invalid uvw shape: must have 3 or 4 values.')

        # Convert to primitive cell
        uvw = miller.vector_conventional_to_primitive(uvw, self.conventional_setting)

        # Convert to Cartesian wrt rotated system
        cart = self.transform.dot(miller.vector_crystal_to_cartesian(uvw, self.ucell.box))

        # Check that Cartesian vector is in slip plane
        if not np.isclose(cart[self.cutindex], 0.0):
            raise ValueError(f'shift vector {value} not in fault plane {self.hkl}')

        # Save uvw and cart
        self.__a2vect_uvw = value
        self.__a2vect_cart = cart

    @property
    def a1vect_cart(self) -> np.ndarray:
        """numpy.ndarray : One of the two shift vectors in Cartesian relative to system."""
        return self.__a1vect_cart
    
    @property
    def a2vect_cart(self) -> np.ndarray:
        """numpy.ndarray : One of the two shift vectors in Cartesian relative to system."""
        return self.__a2vect_cart
    
    @property
    def faultpos_cart(self) -> float:
        """float : The Cartesian position of the slip plane."""
        try:
            return self.__faultpos_cart
        except:
            raise AttributeError('system not yet built. Use build_system() or surface().')
        
    @property
    def faultpos_rel(self) -> float:
        """float : The fractional position of the slip plane."""
        try:
            return self.__faultpos_rel
        except:
            raise AttributeError('system not yet built. Use build_system() or surface().')
        
    @property
    def abovefault(self) -> list:
        """list : Indices of all atoms in system above the slip plane."""
        try:
            return self.__abovefault
        except:
            raise AttributeError('system not yet built. Use build_system() or surface().')
        
    @faultpos_cart.setter
    def faultpos_cart(self, value: float):

        # faultpos_rel = (faultpos_cart - origin) / width (for origin, width || cutindex)
        faultpos_rel = ((value - self.system.box.origin[self.cutindex])
                       / self.system.box.vects[self.cutindex, self.cutindex])
        if faultpos_rel < 0.0 or faultpos_rel > 1.0:
            raise ValueError('faultpos is outside system')
        
        self.__faultpos_rel = faultpos_rel
        self.__faultpos_cart = value
        
        # Identify atoms above fault plane position
        self.__abovefault = self.system.atoms.pos[:, self.cutindex] > (self.faultpos_cart)

    @faultpos_rel.setter
    def faultpos_rel(self, value: float):
        
        if value < 0.0 or value > 1.0:
            raise ValueError('faultpos is outside system')
        
        self.__faultpos_rel = value

        # faultpos_cart = origin + faultpos_rel * width (for origin, width || cutindex) 
        self.__faultpos_cart = (self.system.box.origin[self.cutindex] 
                              + self.faultpos_rel 
                              * self.system.box.vects[self.cutindex, self.cutindex])
        
        # Identify atoms above fault plane position
        self.__abovefault = self.system.atoms.pos[:, self.cutindex] > (self.faultpos_cart)

    def surface(self,
                shift: Optional[npt.ArrayLike] = None,
                vacuumwidth: Optional[float] = None,
                minwidth: Optional[float] = None,
                sizemults: Optional[list] = None,
                even: bool = False,
                faultpos_rel: Optional[float] = None,
                faultpos_cart: Optional[float] = None) -> System:
        """
        Generates the free surface atomic system, which is used as the basis for generating
        the stacking fault configuration(s).
        
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
        faultpos_rel : float, optional
            The position to place the slip plane within the system given as a
            relative coordinate along the out-of-plane direction.  faultpos_rel
            and faultpos_cart cannot both be given.  Default value is 0.5 if 
            faultpos_cart is also not given.
        faultpos_cart : float, optional
            The position to place the slip plane within the system given as a
            Cartesian coordinate along the out-of-plane direction.  faultpos_rel
            and faultpos_cart cannot both be given.
        
        Returns
        -------
        atomman.System
            The free surface atomic system.
        """
        
        super().surface(shift=shift, vacuumwidth=vacuumwidth, minwidth=minwidth,
                             sizemults=sizemults, even=even)
        
        # Set Cartesian fault position
        if faultpos_cart is not None:
            if faultpos_rel is not None:
                raise ValueError('faultpos_rel and faultpos_cart cannot both be given')
            self.faultpos_cart = faultpos_cart
        
        # Set relative fault position
        elif faultpos_rel is not None:
            self.faultpos_rel = faultpos_rel
        
        # Set default fault position
        else:
            self.faultpos_rel = 0.5

        return self.system
    
    def fault(self,
              a1: Optional[float] = None,
              a2: Optional[float] = None,
              outofplane: Optional[float] = None,
              faultshift: Optional[npt.ArrayLike] = None,
              minimum_r: Optional[float] = None,
              a1vect_uvw: Optional[npt.ArrayLike] = None,
              a2vect_uvw: Optional[npt.ArrayLike] = None,
              faultpos_cart: Optional[float] = None,
              faultpos_rel: Optional[float] = None) -> System:
        """
        Generates a fault configuration by displacing all atoms above the slip
        plane.

        Parameters
        ----------
        a1 : float, optional
            The fractional coordinate of a1vect to shift by. 
            Default value is 0.0.
        a2 : float, optional
            The fractional coordinate of a2vect to shift by. 
            Default value is 0.0. 
        outofplane : float, optional
            An out-of-plane shift, given in absolute units.
            Default value is 0.0.
        faultshift : array-like object, optional
            The full shifting vector to displace the atoms above the slip
            plane by.  Cannot be given with a1, a2, or outofplane.
        minimum_r : float, optional
            Specifies the minimum allowed interatomic spacing across the slip
            plane.  If any sets of atoms are closer than this value then the
            outofplane shift is increased.  Default value is None, which
            performs no adjustment.
        a1vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.
            Included here for those wishing to override the values set during
            class initialization.
        a2vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.
            Included here for those wishing to override the values set during
            class initialization.
        faultpos_rel : float, optional
            The position to place the slip plane within the system given as a
            relative coordinate along the out-of-plane direction.  Included
            here for those wishing to override the value set when surface()
            was called.  faultpos_rel and faultpos_cart cannot both be given.
        faultpos_cart : float, optional
            The position to place the slip plane within the system given as a
            Cartesian coordinate along the out-of-plane direction.  Included
            here for those wishing to override the value set when surface()
            was called.  faultpos_rel and faultpos_cart cannot both be given.

        Returns
        -------
        atomman.System
            The atomic configuration with stacking fault shift
        """
        # Update uvws and faultpos if given
        if a1vect_uvw is not None:
            self.a1vect_uvw = a1vect_uvw
        if a2vect_uvw is not None:
            self.a2vect_uvw = a2vect_uvw
        if faultpos_cart is not None:
            if faultpos_rel is not None:
                raise ValueError('faultpos_rel and faultpos_cart cannot both be given')
            self.faultpos_cart = faultpos_cart
        elif faultpos_rel is not None:
            self.faultpos_rel = faultpos_rel

        # Define out of plane unit vector 
        ovect = np.zeros(3)
        ovect[self.cutindex] = 1.0
        
        # Identify the two non-cut indices
        inindex = []
        for i in range(3):
            if i != self.cutindex:
                inindex.append(i)
        
        # Calculate faultshift
        if a1 is not None or a2 is not None or outofplane is not None:
            if faultshift is not None:
                raise ValueError('a1, a2, outofplane cannot be given with faultshift')
            if a1 is None:
                a1 = 0.0
            if a2 is None:
                a2 = 0.0
            if outofplane is None:
                outofplane = 0.0
            faultshift = a1 * self.a1vect_cart + a2 * self.a2vect_cart + outofplane * ovect
        
        # Set default faultshift
        elif faultshift is None:
            faultshift = np.array([0.0, 0.0, 0.0])

        # Shift atoms above the fault by faultshift
        sfsystem = deepcopy(self.system)
        sfsystem.atoms.pos[self.abovefault] += faultshift
        sfsystem.wrap()
        
        # Add additional outofplane shift if necessary
        if minimum_r is not None:
            # Get all atoms within minimum_r of fault position
            top_pos = sfsystem.atoms.pos[self.abovefault]
            top_pos = top_pos[top_pos[:, self.cutindex] <= self.faultpos_cart + minimum_r] 
            bot_pos = sfsystem.atoms.pos[~self.abovefault] # 
            bot_pos = bot_pos[bot_pos[:, self.cutindex] >= self.faultpos_cart - minimum_r]
            if top_pos.shape[0] > 0 and bot_pos.shape[0] > 0:
                dmag_min = minimum_r
                dvect_min = None
                
                for i in range(top_pos.shape[0]):
                    dvect = sfsystem.dvect(bot_pos, top_pos[i])
                    if dvect.shape == (3,):
                        dvect = dvect.reshape(1,3)
                    dmag = np.linalg.norm(dvect, axis=1)
                    i = np.argmin(dmag)
                    if dmag[i] < dmag_min:
                        dmag_min = dmag[i]
                        dvect_min = dvect[i]
                
                if dvect_min is not None:
                    
                    new = (minimum_r**2 - dvect_min[inindex[0]]**2 - dvect_min[inindex[1]]**2)**0.5
                    outofplane = new - dvect_min[self.cutindex]
                    faultshift = outofplane * ovect
                    sfsystem.atoms.pos[self.abovefault] += faultshift
                    sfsystem.wrap()
        
        return sfsystem
    
    def iterfaultmap(self,
                     num_a1: Optional[int] = None,
                     num_a2: Optional[int] = None,
                     outofplane: Optional[float] = None,
                     minimum_r: Optional[float] = None,
                     a1vect_uvw: Optional[npt.ArrayLike] = None,
                     a2vect_uvw: Optional[npt.ArrayLike] = None,
                     faultpos_cart: Optional[float] = None,
                     faultpos_rel: Optional[float] = None
                     ) -> Generator[Tuple[float, float, System], None, None]:
        """
        Iterates over generalized stacking fault configurations associated
        with a 2D map of equally spaced a1, a2 coordinates.

        Parameters
        ----------
        num_a1 : int
            The number of a1 values to generate systems for.  
            Default value is 1 (only generate for a1=0.0).
        num_a2 : int
            The number of a2 values to generate systems for.  
            Default value is 1 (only generate for a2=0.0).
        outofplane : float, optional
            An out-of-plane shift, given in absolute units.
            Default value is 0.0.
        minimum_r : float, optional
            Specifies the minimum allowed interatomic spacing across the slip
            plane.  If any sets of atoms are closer than this value then the
            outofplane shift is increased.  Default value is None, which
            performs no adjustment.
        a1vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.
            Included here for those wishing to override the values set during
            class initialization.
        a2vect_uvw : array-like object, optional
            The crystal vector to use for one of the two shifting vectors.
            Included here for those wishing to override the values set during
            class initialization.
        faultpos_rel : float, optional
            The position to place the slip plane within the system given as a
            relative coordinate along the out-of-plane direction.  Included
            here for those wishing to override the value set when surface()
            was called.  faultpos_rel and faultpos_cart cannot both be given.
        faultpos_cart : float, optional
            The position to place the slip plane within the system given as a
            Cartesian coordinate along the out-of-plane direction.  Included
            here for those wishing to override the value set when surface()
            was called.  faultpos_rel and faultpos_cart cannot both be given.
        
        Yields
        ------
        a1 : float
            The a1 fractional coordinate of a1vect. 
        a2 : float
            The a2 fractional coordinate of a2vect. 
        atomman.System
            The fault configuration associated with the a1, a2 shift.
        """
        # Update uvws and faultpos if given
        if a1vect_uvw is not None:
            self.a1vect_uvw = a1vect_uvw
        if a2vect_uvw is not None:
            self.a2vect_uvw = a2vect_uvw
        if faultpos_cart is not None:
            if faultpos_rel is not None:
                raise ValueError('faultpos_rel and faultpos_cart cannot both be given')
            self.faultpos_cart = faultpos_cart
        elif faultpos_rel is not None:
            self.faultpos_rel = faultpos_rel

        if num_a1 is None:
            num_a1 = 1
        if num_a2 is None:
            num_a2 = 1

        # Construct mesh of regular points
        a1s, a2s = np.meshgrid(np.linspace(0, 1, num_a1, endpoint=False),
                               np.linspace(0, 1, num_a2, endpoint=False))

        for a1, a2 in zip(a1s.flat, a2s.flat):
            yield a1, a2, self.fault(a1=a1, a2=a2, outofplane=outofplane,
                                     minimum_r=minimum_r)