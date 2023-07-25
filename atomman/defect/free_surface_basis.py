# coding: utf-8

# Standard Python libraries
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from .. import Box, System
from ..tools import vect_angle, miller
from ..tools.miller import vector_crystal_to_cartesian

def free_surface_basis(hkl: npt.ArrayLike,
                       box: Optional[Box] = None,
                       cutboxvector: str = 'c',
                       maxindex: Optional[int] = None,
                       return_hexagonal: Optional[bool] = None,
                       return_planenormal: bool = False,
                       conventional_setting: Optional[str] = None
                       ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Generates the uvw box vector orientations for a free surface atomic
    system.  In determining the uvw sets, two sets will be in the specified
    hkl plane and one will be out of the plane.  Uses free surface in-plane
    vector determination algorithm by W. Sun and G. Cedar, Surface Science,
    617, 53-59 (2013) to identify two in-plane vectors and the plane normal.
    The shortest in-plane vector is identified, as well as an out-of-plane
    vector close to the plane normal.  The second in-plane vector is then
    selected to be a shortest in-plane vector that is not parallel to the
    first.

    Parameters
    ----------
    hkl : array-like object
        The free surface plane to generate expressed in either 3 indices
        Miller (hkl) format or 4 indices Miller-Bravais (hkil) format.
    box : atomman.Box, optional
        The box object associated with the unit cell. Used to identify the
        best uvw set for the out-of-plane box vector.  Default value uses a
        cubic box.
    cutboxvector : str, optional
        Specifies which of the three box vectors corresponds to the
        out-of-plane vector.  Default value is c.
    maxindex : int, optional
        Max uvw index value to use in identifying the best uvw set for the
        out-of-plane vector.  If not given, will use the largest absolute
        index between the given hkl and the initial in-plane vector guesses.
    return_hexagonal : bool, optional
        Flag for indicating if the returned vectors are expressed in Miller
        [uvw] format (False) or Miller-Bravais [uvtw] format (True).  The
        Miller-Bravais format is only allowed if box is in the standard
        hexagonal setting: a=b!=c, alpha=beta=90, gamma=120.  Default value is
        False if hkl is given in the 3 indices Miller (hkl) format and True if
        it is given in the 4 indices Miller-Bravais (hkil) format.
    return_planenormal : bool, optional
        If True, the computed Cartesian plane normal will also be returned.
        Default value is False.
    conventional_setting : str, optional
        Allows for rotations of a primitive unit cell box to be determined from
        (hkl) indices specified relative to a conventional unit cell.  Allowed
        settings: 'p' for primitive (no conversion), 'f' for face-centered,
        'i' for body-centered, and 'a', 'b', or 'c' for side-centered.  Default
        behavior is to perform no conversion, i.e. take (hkl) relative to the
        given box.

    Returns
    -------
    uvws : numpy.ndarray
        3x3 array of Miller [uvw] vectors or 3x4 array of Miller-Bravais [uvtw] vectors to rotate the unit cell for a free surface configuration.
    planenormal : numpy.ndarray
        The Cartesian plane normal vector.  Only returned if return_planenormal is True.

    Raises
    ------
    ValueError
        If invalid hkl indices values are given.
    AssertionError
        If the search fails to find any of the three [uvw] rotation vectors.
    """
    # Set default box to be cubic
    if box is None:
        box = Box()

    # Check hkl values
    hkl = np.asarray(hkl)

    # Convert hkil to hkl
    if hkl.shape == (4,):
        if box.ishexagonal():
            hkl = miller.plane4to3(hkl)
            if return_hexagonal is None:
                return_hexagonal = True
        else:
            raise ValueError('Miller-Bravais indices given with non-hexagonal box')
    elif hkl.shape == (3,):
        if return_hexagonal is None:
            return_hexagonal = False
        elif return_hexagonal and not box.ishexagonal():
            raise ValueError('cannot return Miller-Bravais indices for non-hexagonal box')
    else:
        raise ValueError('Invalid hkl indices: must be 3 values or 4')

    if np.allclose(hkl, np.asarray(hkl, dtype=int)):
        hkl = np.asarray(hkl, dtype=int)
    else:
        raise ValueError('hkl indices must be integers')

    # Build conventional box if conventional_setting is given
    if conventional_setting is not None:
        primitive_box = box
        p2c_uvws = miller.vector_conventional_to_primitive(np.identity(3),
                                                           setting=conventional_setting)

        conventional_box = System(box=primitive_box).rotate(p2c_uvws).box
        box = conventional_box

    # Find two in-plane box vectors
    if hkl[0] != 0:
        if hkl[1] != 0:
            if hkl[2] != 0:
                # hkl solution
                m = np.lcm.reduce([hkl[0], hkl[1], hkl[2]]) # pylint: disable=no-member
                s = np.sign(hkl[0] * hkl[1] * hkl[2])
                a_uvw = np.array([-m / hkl[0], m / hkl[1], 0], dtype=int)
                b_uvw = np.array([-m / hkl[0], 0, m / hkl[2]], dtype=int)
            else:
                # hk0 solution
                m = np.lcm(hkl[0], hkl[1]) # pylint: disable=no-member
                s = np.sign(hkl[0] * hkl[1])
                a_uvw = np.array([-m / hkl[0], m / hkl[1], 0], dtype=int)
                b_uvw = np.array([0, 0, 1], dtype=int)
        else:
            if hkl[2] != 0:
                # h0l solution
                m = np.lcm(hkl[0], hkl[2]) # pylint: disable=no-member
                s = np.sign(hkl[0] * hkl[2])
                a_uvw = np.array([m / hkl[0], 0, -m / hkl[2]], dtype=int)
                b_uvw = np.array([0, 1, 0], dtype=int)
            else:
                # h00 solution
                m = 1
                s = np.sign(hkl[0])
                a_uvw = np.array([0, 1, 0], dtype=int)
                b_uvw = np.array([0, 0, 1], dtype=int)
    elif hkl[1] != 0:
        if hkl[2] != 0:
            # 0kl solution
            m = np.lcm(hkl[1], hkl[2]) # pylint: disable=no-member
            s = np.sign(hkl[1] * hkl[2])
            a_uvw = np.array([0, -m / hkl[1], m / hkl[2]], dtype=int)
            b_uvw = np.array([1, 0, 0], dtype=int)
        else:
            # 0k0 solution
            m = 1
            s = np.sign(hkl[1])
            a_uvw = np.array([0, 0, 1], dtype=int)
            b_uvw = np.array([1, 0, 0], dtype=int)
    elif hkl[2] != 0:
        # 00l solution
        m = 1
        s = np.sign(hkl[2])
        a_uvw = np.array([1, 0, 0], dtype=int)
        b_uvw = np.array([0, 1, 0], dtype=int)
    else:
        raise ValueError('hkl cannot be all zeros')

    # Convert in-plane box vectors to primitive cell vectors if needed
    if conventional_setting is not None:
        a_uvw = miller.vector_conventional_to_primitive(a_uvw, setting=conventional_setting)
        b_uvw = miller.vector_conventional_to_primitive(b_uvw, setting=conventional_setting)
        box = primitive_box

    # Set default n if needed
    if maxindex is None:
        maxindex = int(np.max([np.abs(a_uvw), np.abs(b_uvw), np.abs(hkl)]))

    # Compute Cartesian plane normal
    planenormal = s * np.cross(vector_crystal_to_cartesian(a_uvw, box),
                               vector_crystal_to_cartesian(b_uvw, box))

    # Build gen_vector iterator for testing vectors
    def gen_vector(n):
        for kk in range(0, n+1):
            for sk in [1, -1]:
                k = sk * kk
                for jj in range(0, n+1):
                    for sj in [1, -1]:
                        j = sj * jj
                        for ii in range(0, n+1):
                            for si in [1, -1]:
                                i = si * ii
                                if i==0 and j==0 and k==0:
                                    continue
                                yield np.array([i, j, k], dtype=int)

    # First search
    a_mag = np.linalg.norm(vector_crystal_to_cartesian([maxindex, maxindex, maxindex], box))
    c_angle = 90
    a_uvw = None
    c_uvw = None
    for uvw in gen_vector(maxindex):
        cart = vector_crystal_to_cartesian(uvw, box)
        mag = np.linalg.norm(cart)
        angle = vect_angle(cart, planenormal)

        # Find shortest vector in the plane
        if np.isclose(np.dot(cart, planenormal), 0.0):
            if mag < a_mag:
                a_uvw = uvw
                a_mag = mag

        # Find vector closest to plane normal
        elif angle < c_angle:
            c_angle = angle
            c_uvw = uvw

    assert a_uvw is not None, 'Failed to find first vector in slip plane'
    assert c_uvw is not None, 'Failed to find vector near slip plane normal'

    # Reduce c_uvw if possible
    c_uvw = c_uvw / np.gcd.reduce(np.asarray(c_uvw, dtype=int)) # pylint: disable=no-member

    # Second search
    a_cart = vector_crystal_to_cartesian(a_uvw, box)
    b_mag = np.linalg.norm(vector_crystal_to_cartesian([maxindex, maxindex, maxindex], box))
    b_uvw = None
    min_angle = 180.0
    for uvw in gen_vector(maxindex):
        cart = vector_crystal_to_cartesian(uvw, box)
        angle = vect_angle(a_cart, cart)

        # Check that vector is in plane and not parallel to a_uvw
        if np.isclose(np.dot(cart, planenormal), 0.0) and not np.isclose(angle, 0.0) and not np.isclose(angle, 180.0):

            # Check if right-handed
            if np.dot(np.cross(a_cart, cart), planenormal) > 0:

                # Find b_uvw with smallest magnitude and smallest angle
                mag = np.linalg.norm(cart)
                if (np.isclose(mag, b_mag) and angle < min_angle) or mag < b_mag:
                    b_uvw = uvw
                    b_mag = mag
                    min_angle = angle

    assert b_uvw is not None, 'Failed to find second vector in slip plane'

    # Orient the uvw sets based on cutboxvector
    if cutboxvector == 'c':
        uvws = np.array([a_uvw, b_uvw, c_uvw])
    elif cutboxvector == 'b':
        uvws = np.array([b_uvw, c_uvw, a_uvw])
    elif cutboxvector == 'a':
        uvws = np.array([c_uvw, a_uvw, b_uvw])

    if return_hexagonal:
        uvws = miller.vector3to4(uvws)

    if return_planenormal:
        return uvws, planenormal
    else:
        return uvws
