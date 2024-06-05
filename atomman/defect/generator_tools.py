from typing import Tuple

import numpy as np
import numpy.typing as npt

# Local imports
from .. import Box
from ..tools import vect_angle, miller

def planenormal_uvw(box: Box,
                    planenormal: npt.ArrayLike,
                    maxindex: int = 10) -> np.ndarray:
    """
    Finds a Miller crystal vector close to the plane normal.  This algorithm
    first identifies all Miller vectors in the search range that have the
    smallest angle with respect to the plane normal, then filters based on
    the vector magnitudes.
    
    Parameters
    ----------
    box : atomman.Box
        The unit cell box to use as the reference for the Miller crystal
        vectors being searched.
    planenormal : array-like object
        The plane normal expressed as a Cartesian vector in the same system
        as the vectors of box.
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
    # Get unique uvws in the search range
    uvws = miller.all_indices(maxindex=maxindex, reduce=True)
    carts = box.vector_crystal_to_cartesian(uvws)
            
    # Find angles between plane normal and the vectors 
    angles = vect_angle(carts, planenormal)

    # Pick only the vectors with the smallest angle
    at_min_angle = np.isclose(angles, angles.min(), rtol=0, atol=1e-6)
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

def short_inplane_uvw(box: Box,
                      planenormal: npt.ArrayLike,
                      maxindex: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Finds one of the shortest Miller crystal vectors that is in the plane
    defined by a Cartesian plane normal.

    Parameters
    ----------
    box : atomman.Box
        The unit cell box to use as the reference for the Miller crystal
        vectors being searched.
    planenormal : array-like object
        The plane normal expressed as a Cartesian vector in the same system
        as the vectors of box.
    maxindex : int, optional
        The maximum absolute vector index to use for the u v w values:
        u, v, and w can all independently vary from -maxindex to maxindex.
        Default value is 10.
    
    Returns
    -------
    uvw : numpy.ndarray
        A Miller [uvw] crystal vector from the search range that is perpendicular
        to the planenormal and has the shortest magnitude.
    cart : numpy.ndarray
        The Cartesian vector associated with the returned uvw.
    """
    # Get all uvws in the search range
    uvws = miller.all_indices(maxindex=maxindex, reduce=True)
    carts = miller.vector_crystal_to_cartesian(uvws, box)
    
    # Find all that are in-plane (perpendicular to plane normal)
    is_in_plane = np.isclose(np.einsum('...i,...i', planenormal, carts), 0.0)
    uvws = uvws[is_in_plane]
    carts = carts[is_in_plane]
    
    # Find the first shortest vector
    mags = np.linalg.norm(carts, axis=1)
    is_min_mag = np.isclose(mags, mags.min())
    uvw = uvws[is_min_mag][0]
    cart = carts[is_min_mag][0]
    
    return uvw, cart

def second_inplane_uvw(box: Box,
                       planenormal: npt.ArrayLike,
                       uvw1: npt.ArrayLike,
                       maxindex: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Given one in-plane vector, find a second that is near perpendicular and
    results in the smallest 2D cell size.

    Parameters
    ----------
    box : atomman.Box
        The unit cell box to use as the reference for the Miller crystal
        vectors being searched.
    planenormal : array-like object
        The plane normal expressed as a Cartesian vector in the same system
        as the vectors of box.
    uvw1 : array-like object
        A Miller [uvw] crystal vector in the plane.
    maxindex : int, optional
        The maximum absolute vector index to use for the u v w values:
        u, v, and w can all independently vary from -maxindex to maxindex.
        Default value is 10.
    
    Returns
    -------
    uvw : numpy.NDArray
        A Miller [uvw] crystal vector from the search range that forms a
        primitive in-plane 2D cell with uvw1 and is near perpendicular to uvw1.
    cart : numpy.ndarray
        The Cartesian vector associated with the returned uvw.
    """
    uvw1 = np.asarray(uvw1)
    cart1 =  miller.vector_crystal_to_cartesian(uvw1, box)
    
    # Get all uvws in the search range
    uvws = miller.all_indices(maxindex=maxindex, reduce=True)
    carts = miller.vector_crystal_to_cartesian(uvws, box)
    
    # Find all that are in-plane (perpendicular to plane normal)
    is_in_plane = np.isclose(np.einsum('...i,...i', planenormal, carts), 0.0)
    uvws = uvws[is_in_plane]
    carts = carts[is_in_plane]

    # Find all second vectors that give the smallest cell area
    cell_areas = np.linalg.norm(np.cross(cart1, carts), axis=1)
    is_min_area = np.isclose(cell_areas, cell_areas[~np.isclose(cell_areas, 0.0)].min())
    uvws = uvws[is_min_area]
    carts = carts[is_min_area]
    #mags = mags[is_min_area]
    
    # Find angles between the in-plane vectors
    angles = vect_angle(cart1, carts)
    
    # Correct angles for right-handed orientations for uniqueness
    lh = np.dot(np.cross(cart1, carts), planenormal) < 0
    angles[lh] = 360 - angles[lh]
    
    # Identify angles closest to 90
    diff_angle = np.abs(angles - 90)
    is_min_diff_angle = np.isclose(diff_angle, diff_angle.min())
    
    # Select angle < 90 for symmetric cases
    best_angle = angles[is_min_diff_angle].min()
    is_best_angle = np.isclose(best_angle, angles)
    
    uvw2 = uvws[is_best_angle][0]
    cart2 = carts[is_best_angle][0]
    
    return uvw2, cart2
    