# coding: utf-8
# Standard Python imports
from typing import Optional, Tuple
from copy import deepcopy
from math import ceil

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from scipy.spatial.transform import Rotation

# Local imports
from .. import System
from ..tools import vect_angle, miller
from .generator_tools import second_inplane_uvw, planenormal_uvw

class TiltGrainBoundaryHelper():
    """
    Supporting class that collects various calculation tools for tilt grain
    boundaries.
    """
    
    def __init__(self,
                 ucell: System,
                 axis_uvw: npt.ArrayLike,
                 conventional_setting: str = 'p',
                 ref_uvw: Optional[npt.ArrayLike] = None,
                 ref_hkl: Optional[npt.ArrayLike] = None,
                 tol: float = 1e-8):
        """
        Class initializer.  Takes inputs that allow for the definition of the
        tilt grain boundary family (unit cell and tilt axis).
        
        Parameters
        ----------
        ucell : am.System
            The unit cell to use as the reference for the grains.
        axis_uvw : array-like object
            The Miller [uvw] or Miller-Bravais [uvtw] crystal vector of ucell
            that serves as the tilt axis.  The tilt axis is a common vector in
            the grain boundary plane that is shared by both grains.
        conventional_setting : str, optional
            Specifies which conventional lattice setting that ucell is in.
            The default value of 'p' takes ucell to be primitive, in which case
            uvw input values are limited to integers.  The correct setting is
            needed in order to use non-integer lattice vectors.
        ref_hkl : array-like object, optional
            A Miller (hkl) or Miller-Bravais (hkil) plane to use as the
            reference 0 degree tilt plane. This is used by the symmetric tilt
            finder to properly and consistently find the correct second grain
            orientation.  Alternatively, the 0 degree tilt plane can be
            specified using ref_uvw.  If neither ref_hkl nor ref_uvw are given
            then the second grain's boundary plane is selected based on
            producing the smallest misorientation angle which often, but not
            always, corresponds to a symmetric tilt boundary.
        ref_uvw : array-like object, optional
            A Miller [uvw] or Miller-Bravais [uvtw] vector in the grain boundary
            plane that along with axis_uvw identifies the reference 0 degree
            tilt plane.  This is used by the symmetric tilt
            finder to properly and consistently find the correct second grain
            orientation.  Alternatively, the 0 degree tilt plane can be
            specified using ref_hkl.  If neither ref_hkl nor ref_uvw are given
            then the second grain's boundary plane is selected based on
            producing the smallest misorientation angle which often, but not
            always, corresponds to a symmetric tilt boundary.
        tol : float, optional
            A floating point tolerance used in finding the compatible primitive
            unit cell.  Available as a parameter in case issues are encountered.
            Default value is 1e-8.
        """
        
        # Initial handling of ucell
        if not isinstance(ucell, System):
            raise TypeError('ucell must be an atomman.System object')
        
        # Create compatible primitive unit cell
        ucell_prim, transform_c_to_p = ucell.dump('conventional_to_primitive',
                                                  setting=conventional_setting,
                                                  return_transform=True, atol=tol)
        
        # Save unit cells and transformation properties
        self.__ucell = ucell
        self.__ucell_prim = ucell_prim
        self.__transform_c_to_p = Rotation.from_matrix(transform_c_to_p)
        self.__conventional_setting = conventional_setting
        
        # Initial handling of axis_uvw
        axis_uvw = np.asarray(axis_uvw)
        axis_uvw_prim, axis_cart_prim = self.primitive_vectors(axis_uvw)
        
        # Save vectors
        self.__axis_uvw = axis_uvw
        self.__axis_uvw_prim = axis_uvw_prim
        self.__axis_cart_prim = axis_cart_prim

        # Convert ref_hkl into a Cartesian plane normal wrt ucell_prim
        if ref_hkl is not None:
            if ref_uvw is not None:
                raise ValueError('ref_hkl and ref_uvw cannot be given together')
            ref_hkl = np.asarray(ref_hkl)
            ref_plane_norm = self.ucell.box.plane_crystal_to_cartesian(ref_hkl)
            ref_plane_norm_prim = self.transform_c_to_p.apply(ref_plane_norm)

        # Convert ref_uvw to Cartesian wrt ucell_prim and cross with axis
        elif ref_uvw is not None:
            ref_uvw = np.asarray(ref_uvw)
            ref_uvw_prim, ref_cart_prim = self.primitive_vectors(ref_uvw)
            ref_plane_norm_prim = np.cross(self.axis_cart_prim, ref_cart_prim)

        # Set default None value otherwise
        else:
            ref_plane_norm_prim = None

        # Save reference plane normal
        self.__ref_plane_norm_prim = ref_plane_norm_prim

    @property
    def ucell(self) -> System:
        """atomman.System: The conventional crystal unit cell"""
        return self.__ucell
    
    @property
    def ucell_prim(self) -> System:
        """atomman.System: The primitive crystal unit cell"""
        return self.__ucell_prim
    
    @property
    def conventional_setting(self) -> str:
        """str: The conventional lattice setting associated with ucell"""
        return self.__conventional_setting

    @property
    def transform_c_to_p(self) -> Rotation:
        """
        scipy.spatial.transform.Rotation: The rotation transformation for changing from
        ucell's Cartesian system to ucell_prim's.        
        """
        return self.__transform_c_to_p

    @property
    def axis_uvw(self) -> np.ndarray:
        """numpy.NDArray: The tilt axis as a Miller or Miller-Bravais crystal vector of ucell""" 
        return self.__axis_uvw
    
    @property
    def axis_uvw_prim(self) -> np.ndarray:
        """numpy.NDArray: The tilt axis as a Miller crystal vector of ucell_prim""" 
        return self.__axis_uvw_prim

    @property
    def axis_cart_prim(self) -> np.ndarray:
        """numpy.NDArray: The tilt axis as a Cartesian vector wrt ucell_prim's orientation""" 
        return self.__axis_cart_prim
    
    @property
    def ref_plane_norm_prim_prim(self) -> Optional[np.ndarray]:
        """numpy.NDArray or None: The plane normal wrt ucell_prim for the 0 degree reference plane""" 
        return self.__ref_plane_norm_prim

    def symmetric_uvws(self,
                       plane1_hkl: Optional[npt.ArrayLike] = None,
                       in1_uvw: Optional[npt.ArrayLike] = None,
                       cutboxvector: str = 'c',
                       maxindex: int = 20,
                       ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Given an in-plane vector for grain #1 (not parallel to axis),
        find the uvw crystal vector rotation sets for the two grains that
        generates a symmetric tilt boundary based in that in-plane vector and
        the tilt axis.

        Parameters
        ----------
        plane1_hkl : array-like object, optional
            The Miller (hkl) or Miller-Bravais (hkil) crystal plane of ucell
            that is the boundary plane for the first grain.  Either plane1_hkl
            or in1_uvw must be given to identify the tilt boundary.
        in1_uvw : array-like object, optional
            The second in-plane Miller [uvw] or Miller-Bravais [uvtw] crystal
            vector of ucell that along with axis_uvw defines the boundary plane
            for the first grain.  Either plane1_hkl or in1_uvw must be given to
            identify the tilt boundary.
        cutboxvector : str, optional
            Sets the alignment of the uvws to the final configuration's box
            vectors by specifying which of the box vectors is out-of-plane
            with respect to the grain boundary plane. Default value is 'c',
            which aligns the tilt axis with a, and in1 (and in2) with b.
        maxindex : int, optional
            The maximum absolute vector index to search over for the u,v,w
            values. u, v, and w can all independently vary from
            -maxindex to maxindex.  Default value is 20.

        Returns
        -------
        uvws1 : numpy.ndarray
            The three Miller crystal vectors that specify how the unit cell
            is to be rotated in order to form grain #1.
        uvws2 : numpy.ndarray
            The three Miller crystal vectors that specify how the unit cell
            is to be rotated in order to form grain #2.
        misorientation : float
            The misorientation angle in degrees for the associated symmetric
            tilt boundary.
        """

        # Convert in1_uvw to the primitive cell and find plane normal
        if in1_uvw is not None:
            if plane1_hkl is not None:
                raise ValueError('in1_uvw and plane1_hkl cannot both be given')
            in1_uvw = np.asarray(in1_uvw, dtype=float)
            in1_uvw_prim, in1_cart_prim = self.primitive_vectors(in1_uvw)
            plane1_norm_prim = np.cross(self.axis_cart_prim, in1_cart_prim)

        # Find plane normal and second in-plane vector using plane1_hkl
        elif plane1_hkl is not None:
            plane1_norm = self.ucell.box.plane_crystal_to_cartesian(plane1_hkl)
            plane1_norm_prim = self.transform_c_to_p.apply(plane1_norm)
            in1_uvw_prim, in1_cart_prim = second_inplane_uvw(self.ucell_prim.box,
                                                             plane1_norm_prim,
                                                             self.axis_uvw_prim,
                                                             maxindex=maxindex)

        else:
            raise ValueError('in1_uvw or plane1_hkl must be given')

        # Find the out of plane uvw closest to the plane normal for grain #1
        pn1_uvw_prim = planenormal_uvw(self.ucell_prim.box, plane1_norm_prim, maxindex=maxindex)

        # Compute |in1| and the angle between in1 and axis
        in1_mag = np.linalg.norm(in1_cart_prim)
        axis_in1_angle = vect_angle(in1_cart_prim, self.axis_cart_prim)

        # Generate a list of all guess in2 values
        in2_uvw_prims = miller.all_indices(maxindex)
        in2_cart_prims = self.ucell_prim.box.vector_crystal_to_cartesian(in2_uvw_prims)
        
        # Remove uvws parallel to in1
        parallel = np.all(np.isclose(np.cross(in2_cart_prims, in1_cart_prim), 0.0, atol=1e-5, rtol=0.0), axis=1)
        in2_uvw_prims = in2_uvw_prims[~parallel]
        in2_cart_prims = in2_cart_prims[~parallel]

        # Limit to only those at the same angle wrt axis
        axis_in2_angles = vect_angle(in2_cart_prims, self.axis_cart_prim)
        same_angle = np.isclose(axis_in2_angles, axis_in1_angle)
        in2_uvw_prims = in2_uvw_prims[same_angle]
        in2_cart_prims = in2_cart_prims[same_angle]
        
        # Limit to only values with same mag
        in2_mags = np.linalg.norm(in2_cart_prims, axis=1)
        same_mag = np.isclose(in1_mag, in2_mags)
        in2_uvw_prims = in2_uvw_prims[same_mag]
        in2_cart_prims = in2_cart_prims[same_mag]

        # Compute the misorientation angles between plane norms
        plane2_norm_prims = np.cross(self.axis_cart_prim, in2_cart_prims)
        misorientations = vect_angle(plane2_norm_prims, plane1_norm_prim)

        # Identify symmetric tilt using ref plane
        if self.ref_plane_norm_prim_prim is not None:
            half_misorientation1 = vect_angle(plane1_norm_prim, self.ref_plane_norm_prim_prim)
            
            # Unique case for 0 degree tilt plane
            if np.isclose(half_misorientation1, 0.0):
                in2_uvw_prim = in1_uvw_prim
                in2_cart_prim = in1_cart_prim

            # Find in2 based on full and half misorientations
            else:
                # Check for correct full misorientation between the two gb planes
                misorientation = 2 * half_misorientation1
                right_mis = np.isclose(misorientations, misorientation)

                # Check the half misorientations wrt the reference plane
                half_misorientations2 = vect_angle(plane2_norm_prims, self.ref_plane_norm_prim_prim)
                same_half_mis = np.isclose(half_misorientations2, half_misorientation1)
                in2_uvw_prim = in2_uvw_prims[(right_mis) & (same_half_mis)][0]
                in2_cart_prim = in2_cart_prims[(right_mis) & (same_half_mis)][0]

        # Guess for a symmetric tilt using smallest misorientation
        else:
            
            # Apply right-handed fix for angles > 180 for uniqueness
            lh = np.dot(np.cross(plane2_norm_prims, plane1_norm_prim), self.axis_cart_prim) < 0
            misorientations[lh] = 360 - misorientations[lh]

            # Find in2 with smallest misorientation
            misorientation = misorientations.min()
            is_min = misorientations == misorientation
            in2_uvw_prim = in2_uvw_prims[is_min][0]
            in2_cart_prim = in2_cart_prims[is_min][0]

        # Find the out of plane uvw closest to the plane normal for grain #2
        plane2_norm_prim = np.cross(self.axis_cart_prim, in2_cart_prim)
        pn2_uvw_prim = planenormal_uvw(self.ucell_prim.box, plane2_norm_prim, maxindex=maxindex)

        # Build uvws based on cutboxvector
        if cutboxvector == 'a':
            uvws1_prim = np.array([pn1_uvw_prim, self.axis_uvw_prim, in1_uvw_prim])
            uvws2_prim = np.array([pn2_uvw_prim, self.axis_uvw, in2_uvw_prim])
        elif cutboxvector == 'b':
            uvws1_prim = np.array([in1_uvw_prim, pn1_uvw_prim, self.axis_uvw])
            uvws2_prim = np.array([in2_uvw_prim, pn2_uvw_prim, self.axis_uvw])
        elif cutboxvector == 'c':
            uvws1_prim = np.array([self.axis_uvw, in1_uvw_prim, pn1_uvw_prim])
            uvws2_prim = np.array([self.axis_uvw, in2_uvw_prim, pn2_uvw_prim])
        else:
            raise ValueError("cutboxvector values limited to 'a', 'b', or 'c'")

        # Convert uvws to conventional ucell
        uvws1 = miller.vector_primitive_to_conventional(uvws1_prim)
        uvws2 = miller.vector_primitive_to_conventional(uvws2_prim)

        return uvws1, uvws2

    

    def asymmetric_in2s(self,
                        in1_uvw: npt.ArrayLike,
                        #ref_uvw: Optional[npt.ArrayLike] = None,
                        maxindex: Optional[int] = None,):
        """
        Find all in2 Miller vectors within a maxindex range that form a
        compatible tilt boundary with a given in1 Miller vector.

        Parameters
        ----------
        in1_uvw : array-like object
            The first grain's in-plane vector given as a Miller or
            Miller-Bravais crystal vector.
        maxindex : int, optional
            The maximum absolute vector index to use for the u v w values:
            u, v, and w can all independently vary from -maxindex to maxindex.
            Default value is the greater of 30 or 4 times the max index within
            in1_uvw.
        """
        # Get axis vector and ucell's box
        box = self.ucell.box
        axis_cart = self.axis_cart
        
        # Initial handling of in1_uvw
        in1_uvw = np.asarray(in1_uvw)
        if in1_uvw.shape == (4,):
            in1_uvw = miller.vector4to3(in1_uvw)
        elif in1_uvw.shape != (3,):
            raise ValueError('invalid dimensions for in1_uvw: must be (3,) or (4,)')
        in1_cart = miller.vector_crystal_to_cartesian(in1_uvw, box)
        
        # Compute |in1|, angle between in1 and axis, and the in1-axis plane normal
        in1_mag = np.linalg.norm(in1_cart)
        axis_in1_angle = vect_angle(in1_cart, axis_cart)
        axis_in1_plane_norm = np.cross(in1_cart, axis_cart)

        # Set default maxindex if nothing is given
        if maxindex is None:
            maxindex = 4 * np.max(np.abs(in1_uvw))
            if maxindex < 30:
                maxindex = 30
        
        # Generate a list of all guess in2 values
        in2_uvws = miller.all_indices(maxindex)
        in2_carts = box.vector_crystal_to_cartesian(in2_uvws)
            
        # Exclude in2 || in1
        parallel_to_in1 = np.all(np.isclose(np.cross(in2_carts, in1_cart), 0.0, atol=1e-5, rtol=0.0), axis=1)
        in2_uvws = in2_uvws[~parallel_to_in1]
        in2_carts = in2_carts[~parallel_to_in1]
        
        # Only use in2 with same angle wrt axis as in1
        axis_in2_angles = vect_angle(in2_carts, axis_cart)
        same_angle = np.isclose(axis_in2_angles, axis_in1_angle)
        in2_uvws = in2_uvws[same_angle]
        in2_carts = in2_carts[same_angle]

        # Only use in2 with compatible mags with in1
        in2_mags = np.linalg.norm(in2_carts, axis=1)
        compatible_mag = np.isclose(in1_mag, in2_mags)

        maxmult = int(ceil(maxindex / np.max(np.abs(in1_uvw)))) + 1
        for mult in range(2, maxmult):
            compatible_mag = compatible_mag | np.isclose(mult * in1_mag, in2_mags)

        in2_uvws = in2_uvws[compatible_mag]
        in2_carts = in2_carts[compatible_mag]
        in2_mags = in2_mags[compatible_mag]

        # Compute the misorientation angles between plane norms
        axis_in2_plane_norms = np.cross(in2_carts, axis_cart)
        misorientations = vect_angle(axis_in2_plane_norms, axis_in1_plane_norm)

        # Apply right-handed fix for angles > 180 for uniqueness
        lh = np.dot(np.cross(axis_in2_plane_norms, axis_in1_plane_norm), axis_cart) < 0
        misorientations[lh] = 360 - misorientations[lh]
        
        # Find only the shortest uvws for each unique misorientation
        small_indices = []
        for misorientation in np.unique(misorientations):
            same_mis = np.isclose(misorientation, misorientations)
            min_mag = in2_mags[same_mis].min()
            small_index = np.where(same_mis & np.isclose(in2_mags, min_mag))[0][0]
            if small_index not in small_indices:
                small_indices.append(small_index)
        in2_uvws = in2_uvws[small_indices]
        in2_carts = in2_carts[small_indices]
        misorientations = misorientations[small_indices]

        # Reduce size if possible
        for i in range(in2_uvws.shape[0]):
            gcd = np.gcd.reduce(in2_uvws[i])
            in2_uvws[i] = in2_uvws[i] / gcd
            in2_carts[i] = in2_carts[i] / gcd

        return in2_uvws, in2_carts, misorientations


    
    def primitive_vectors(self, uvw: np.ndarray):
        """
        Utility function for crystal vector input values to convert them from
        the conventional setting to the primitive, and to compute the
        associated Cartesian values

        Parameters
        ----------
        uvw : numpy.ndarray
            A Miller [uvw] or Miller-Bravais [utww] crystal vector as defined
            relative to the conventional unit cell used to initialize the
            helper object.

        Returns
        -------
        uvw_prim : numpy.ndarray
            The input vector transformed into a Miller [uvw] crystal vector
            relative to the primitive unit cell.
        cart_prim : numpy.ndarray
            The input vector transformed into a Cartesian vector relative to
            the primitive unit cell's Cartesian system.
        """
        # Change [uvtw] into a [uvw]
        if uvw.shape == (4,):
            uvw = miller.vector4to3(uvw)
        elif uvw.shape != (3,):
            raise ValueError('invalid dimensions for Miller(-Bravais) vector: must be (3,) or (4,)')
        
        # Convert to the primitive cell
        uvw_prim = miller.vector_conventional_to_primitive(uvw, self.conventional_setting)

        # Convert to Cartesian
        cart_prim = self.ucell_prim.box.vector_crystal_to_cartesian(uvw_prim)

        return uvw_prim, cart_prim