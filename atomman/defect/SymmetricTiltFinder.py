# coding: utf-8
# Standard Python imports
from typing import Optional, Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

import pandas as pd

# Local imports
from .. import System
from ..tools import vect_angle

class SymmetricTiltFinder():
    """Class for identifying symmetric tilt boundary orientations"""
    
    def __init__(self,
                 ucell: System,
                 hkl: npt.ArrayLike,
                 righthand: bool = True,
                 conventional_setting: str = 'p'):
        """
        Initializes an object of the class to identify parameters associated with 
        symmetric tilt grain boundaries for a given unit cell and crystal plane.
        
        Parameters
        ----------
        ucell : am.System
            The unit cell of the crystal that the tilt boundary search is being performed on.
        hkl : array-like
            The Miller hkl or Miller-Bravais hkil crystal plane of ucell that the tilt boundaries
            will be identified for.
        righthand : bool, optional
            If True, a right-handed convention will be used in determining pairs of in-plane vectors
            and the angles between them.  If False, the left-handed pairs will be returned instead.
            If you get different tilt angles than you expect for a given vector, try changing this
            setting.
        conventional_setting: str, optional
            Indicates the space lattice setting of the given unit cell, i.e.
            'p' for primitive, 'i' for body-centered, 'f' for face-centered,
            'a', 'b', or 'c' for side-centered and 't1', or 't2' for trigonal
            in a hexagonal setting.  Setting this with the appropriate
            conventional unit cell allows for identifying lattice vectors that
            are not integers with respect to the conventional unit cell.  This
            also creates the rotated cell from a compatible primitive cell,
            thereby the final dislocation configurations can be smaller than
            possible solely from the conventional unit cell.
        """
        
        # Set unit cell
        assert isinstance(ucell, System)
        ucell_primitive, transform_conv_prim = ucell.dump('conventional_to_primitive',
                                                          setting=conventional_setting,
                                                          return_transform=True)
        self.__ucell = ucell
        self.__conventional_setting = conventional_setting
        self.__ucell_primitive = ucell_primitive
        self.__transform_conv_prim = transform_conv_prim
        
        # Set hkl plane
        hkl = np.asarray(hkl)
        self.__planenormal = ucell.box.plane_crystal_to_cartesian(hkl)
        self.__hkl = hkl 
        
        # Set righthand flag
        self.righthand = righthand
        
    @property
    def ucell(self) -> System:
        """atomman.System: The crystal unit cell"""
        return self.__ucell
    
    @property
    def conventional_setting(self) -> str:
        """str: The space lattice setting for ucell"""
        return self.__conventional_setting
    
    @property
    def ucell_primitive(self) -> System:
        """atomman.System: The compatible primitive cell associated with ucell"""
        return self.__ucell_primitive
    
    @property
    def transform_conv_prim(self) -> np.ndarray:
        """numpy.NDArray: The transformation matrix associated with ucell to ucell_primitive's Cartesian systems"""
        return self.__transform_conv_prim
    
    @property
    def hkl(self) -> np.ndarray:
        """numpy.NDArray: The Miller or Miller-Bravais crystal plane for the boundary""" 
        return self.__hkl
    
    @property
    def planenormal(self) -> np.array:
        """numpy.NDArray: The Cartesian vector wrt ucell for the plane normal"""
        return self.__planenormal
    
    @property
    def righthand(self) -> bool:
        """bool: Flag indicating if right-handed or left-handed tilt vector pairs are returned"""
        return self.__righthand
    
    @righthand.setter
    def righthand(self, value: bool):
        if not isinstance(True, bool):
            raise TypeError('righthand must be a bool')
        self.__righthand = value
    
    def all_uvws(self,
                maxindex: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds all Miller uvw vectors with u v w indices less than maxindex.

        Parameters
        ----------
        maxindex : int, optional
            The maximum absolute vector index to use for the u v w values:
            u, v, and w can all independently vary from -maxindex to maxindex.
            Default value is 10.

        Returns 
        -------
        uvws : numpy.NDArray
            All [uvw] Miller vectors within the maxindex range.
        carts : numpy.NDArray
            The corresponding Cartesian vectors for all the uvws Miller vectors. 
        """
        # Indices range -maxindex to maxindex
        indices = np.arange(-maxindex, maxindex+1)
        
        # Build grid, combine and transform
        u,v,w = np.meshgrid(indices, indices, indices)
        uvws = np.vstack([u.flat,v.flat,w.flat]).T
        
        # Remove [0,0,0]
        uvws = uvws[~np.all(np.isclose(uvws, 0.0),axis=1)]

        # Convert into Cartesian vectors
        carts = self.ucell.box.vector_crystal_to_cartesian(uvws)

        return uvws, carts

    def planenormal_uvw(self,
                        maxindex: int = 10) -> np.ndarray:
        """
        Finds a Miller crystal vector close to the plane normal.  This algorithm
        first identifies all Miller vectors in the search range that have the
        smallest angle with respect to the plane normal, then filters based on
        the vector magnitudes.
        
        Parameters
        ----------
        maxindex : int, optional
            The maximum absolute vector index to use for the u v w values:
            u, v, and w can all independently vary from -maxindex to maxindex.
            Default value is 10.
        
        Returns
        -------
        uvw : numpy.NDArray
            A Miller [uvw] crystal vector from the search range that has the
            smallest angle with the plane normal.
        """
        # Get all uvws in the search range
        uvws, carts = self.all_uvws(maxindex=maxindex)
                
        # Find angles between plane normal and the vectors 
        angles = vect_angle(carts, self.planenormal)

        # Pick only the vectors with the smallest angle
        at_min_angle = np.isclose(angles, np.min(angles), rtol=0, atol=1e-6)
        uvws = uvws[at_min_angle]
        carts = carts[at_min_angle]

        # Select one of the shortest vectors at min angle
        if uvws.shape[0] == 1:
            uvw = uvws[0]
            
        elif uvws.shape[0] > 0:
            mags = np.linalg.norm(carts, axis=1)
            uvw = uvws[np.isclose(mags, mags.min())][0]
        
        else:
            raise RuntimeError('no angles matching min(angle)!?')
            
        return uvw 

    def inplane_uvws(self,
                     maxindex: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds all Miller uvw vectors with u v w indices less than maxindex
        that are in the grain boundary plane.
        
        Parameters
        ----------
        maxindex : int, optional
            The maximum absolute vector index to use for the u v w values:
            u, v, and w can all independently vary from -maxindex to maxindex.
            Default value is 10.
            
        Returns
        ------
        uvws : numpy.NDArray
            The Miller [uvw] crystal vectors of the in-plane vectors.
        carts : numpy.NDArray
            The corresponding Cartesian vectors.
        """
        uvws, carts = self.all_uvws(maxindex=maxindex)
        
        # Define in-plane filter
        inplane = np.isclose(np.dot(carts, self.planenormal), 0.0)
        
        return uvws[inplane], carts[inplane]
       
    def rotation_uvws(self,
                      uvw1: npt.ArrayLike,
                      maxindex: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Identifies the [uvw] Miller crystal vectors to use in rotating the
        unit cell to generate the appropriate tilt boundary.

        NOTE: Currently, the returned uvws are hard-coded such that the
        out-of-plane direction corresponds to the c box vector. A setting to
        change this is planned to be added.

        Parameters
        ----------
        uvw1 : Array-like object
            The in-plane [uvw] Miller crystal vector that combined with the
            (hkl) plane defines the tilt boundary.
        maxindex : int or None, optional
            The maximum absolute vector index to use for the u v w values:
            u, v, and w can all independently vary from -maxindex to maxindex.
            If None (default) then the largest index in uvw1 and the hkl plane
            definition will be used.

        Returns
        -------
        uvws1 : numpy.NDArray
            The three [uvw] Miller rotation vectors to use on the unit cell to
            generate one grain in the grain boundary configuration.
        uvws2 : numpy.NDArray
            The three [uvw] Miller rotation vectors to use on the unit cell to
            generate the other grain in the grain boundary configuration.
        angle : float
            The tilt boundary angle.
        """
        # Convert and check uvw1 (4 to 3 conversion needed?)
        uvw1 = np.asarray(uvw1)
        cart1 = self.ucell.box.vector_crystal_to_cartesian(uvw1)
        mag1 = np.linalg.norm(cart1)

        # Set default maxindex if nothing is given
        if maxindex is None:
            maxindex = np.max(np.abs(uvw1))

        # Get all in-plane uvws, Cartesian vectors, mags and angles wrt uvw1
        uvws, carts = self.inplane_uvws(maxindex)
        mags = np.linalg.norm(carts, axis=1)
        angles_wrt_1 = vect_angle(carts, cart1)

        # fix angles > 180 based on right/left-hand orientations
        if self.righthand:
            lh = np.dot(np.cross(cart1, carts), self.planenormal) < 0
            angles_wrt_1[lh] = 360 - angles_wrt_1[lh]
        else:
            rh = np.dot(np.cross(cart1, carts), self.planenormal) > 0
            angles_wrt_1[rh] = 360 - angles_wrt_1[rh]

        # Identify all in-plane vectors of the same length as uvw1
        samemag = np.isclose(mags, mag1) & ~np.all(np.isclose(uvws, uvw1), axis=1)
        
        # Find uvw with the shortest angle and same length
        angle = angles_wrt_1[samemag].min()
        uvw2 = uvws[samemag & (angles_wrt_1 == angle)][0]
        cart2 = carts[samemag & (angles_wrt_1 == angle)][0]

        # Compute in-plane angles with respect to uvw2
        angles_wrt_2 = angles_wrt_1 - angle
        angles_wrt_2[angles_wrt_2 < 0] = angles_wrt_2[angles_wrt_2 < 0] + 360

        # Search for secondary matching in-plane vectors for the grains
        unique_mags = np.unique(np.around(mags, decimals=6))
        found = False
        for unique_mag in unique_mags:
            with_mag = np.isclose(mags, unique_mag, rtol=0, atol=1e-6)
            angles_wrt_1_with_mag = angles_wrt_1[with_mag]
            angles_wrt_2_with_mag = angles_wrt_2[with_mag]

            for unique_angle in np.sort(angles_wrt_1_with_mag):
                if np.isclose(unique_angle, 0.0, rtol=0, atol=1e-6):
                    continue
                if np.sum(np.isclose(angles_wrt_2_with_mag, unique_angle, rtol=0, atol=1e-6)) > 0:
                    found = True
                    break
            if found:
                break

        if not found:
            raise ValueError('no compatible secondary in-plane vectors found')
        
        # Pull out matching secondary vectors
        uvw_b_1 = uvws[np.isclose(mags, unique_mag, rtol=0, atol=1e-6) & np.isclose(angles_wrt_1, unique_angle, rtol=0, atol=1e-6)][0]
        uvw_b_2 = uvws[np.isclose(mags, unique_mag, rtol=0, atol=1e-6) & np.isclose(angles_wrt_2, unique_angle, rtol=0, atol=1e-6)][0]
        
        # Get uvw closest to the plane normal
        uvw_c = self.planenormal_uvw(maxindex=maxindex)

        # Construct the two rotation uvw sets
        uvws1 = np.array([uvw1, uvw_b_1, uvw_c])
        uvws2 = np.array([uvw2, uvw_b_2, uvw_c])

        return uvws1, uvws2, angle
        