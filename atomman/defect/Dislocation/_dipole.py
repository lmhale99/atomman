# coding: utf-8
from typing import Optional, Tuple, Union
from copy import deepcopy

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
from ... import System
from . import VolterraDislocation

@staticmethod
def dipole_displacement(dislsol: VolterraDislocation,
                        pos: npt.ArrayLike,
                        x1: Union[float, npt.ArrayLike],
                        x2: Union[float, npt.ArrayLike],
                        mvect: npt.ArrayLike,
                        nvect: npt.ArrayLike,
                        N: int = 5) -> np.ndarray:
    """
    Uses the method of Cai, Bulatov, Chang, Li & Yip, Phil Mag 2003, 83(5), 539-567
    https://doi.org/10.1080/0141861021000051109 to compute the displacement associated
    with a dislocation dipole in a periodic lattice.

    Parameters
    ----------
    dislsol : atomman.defect.VolterraDislocation
        A Volterra dislocation solution object. This provides the dislocation
        orientation (m, n vectors) and the displacement solutions for the two
        individual dislocations of opposite Burgers vectors.
    pos : array-like object
        The coordinates to evaluate the displacement at.
    x1 : float or array-like object
        The coordinate(s) of the first dislocation.  If a float is given, it will be
        multiplied by the dislocation solution's m vector.
    x2 : float or array-like object
        The coordinate(s) of the second dislocation.  If a float is given, it will be
        multiplied by the dislocation solution's m vector.
    mvect : array-like object
        One of the two 2D cell vectors.  This should be parallel to the dislocation
        solution's m vector and a whole periodic lattice vector.
    nvect : array-like object
        One of the two 2D cell vectors.  This should be a whole lattice vector with
        major component along the dislocation solution's n vector. It does not need
        to be parallel to n and can have components in either of the two coordinate
        directions.
    N : int, optional
        Indicates how many image cells are used. A rectangular grid of images is used
        meaning that there will be NxN total cells evaluated, of which NxN-1 will be
        image cells.  Default value is 5.

    Returns
    -------
    displacement : numpy.NDArray
        The associated displacement field evaluated at pos for the dipole configuration.
    """

    # Extract the dislocation solution m,n axes
    m = dislsol.m
    n = dislsol.n

    # Manage input parameters
    pos = np.asarray(pos)
    mvect = np.asarray(mvect)
    nvect = np.asarray(nvect)
    if isinstance(x1, (float, int)):
        x1 = x1 * m
    else:
        x1 = np.asarray(x1)
    if isinstance(x2, (float, int)):
        x2 = x2 * m
    else:
        x2 = np.asarray(x2)

    # Define the non-periodic component of the displacement
    def nonperiodic(pos):
        """
        The dipole displacement field in an infinite medium, i.e. no periodic
        replicas only the two dislocations themselves.
        """
        disp1 = dislsol.displacement(pos - x1)
        disp2 = dislsol.displacement(pos - x2)

        return disp1 - disp2

    # Define the periodic component of the displacement
    def periodic(pos):
        """
        The dipole displacement field including periodic replicas from i,j = -N to N
        in both directions normal to the dislocation line.  NOTE that this includes the
        non-periodic component i=j=0. 
        """
        disp = np.zeros_like(pos)
        for i in range(-N, N+1):
            for j in range(-N, N+1):
                shifted_pos = pos + i * mvect + j * nvect
                disp += nonperiodic(shifted_pos)

        return disp

    # Get m, n components of mvect and nvect
    mn = np.array([mvect, nvect]).dot(np.array([m, n]).T)

    # Get reciprocal of mn
    reciprocal_mn = np.linalg.inv(mn).T

    # Compute the displacement values at the four corners of the 2D [mvect, nvect] cell
    corners = np.array([
        -0.5 * mvect - 0.5 * nvect, # bottom left corner, blc
         0.5 * mvect - 0.5 * nvect, # blc + mvect
        -0.5 * mvect + 0.5 * nvect, # blc + nvect
         0.5 * mvect + 0.5 * nvect]) # blc + mvect + nvect
    disp_corners = periodic(corners)

    # Find total changes in displacement along mvect and nvect
    delta_disp = np.array([disp_corners[1] - disp_corners[0], disp_corners[2] - disp_corners[0]])

    # Define the correction component of the displacement
    def correction(pos):
        """
        This computes the linear displacement correction to apply to the system
        to ensure that it is compatible across the periodic boundaries
        """
        # Transform pos to be relative to mn and multiply by delta disp
        pos_mn = np.vstack([pos.dot(m), pos.dot(n)]).T
        relpos_mn = np.inner(pos_mn, reciprocal_mn)
        disp = relpos_mn.dot(delta_disp)

        return disp

    # Verify the correction against the corners
    correction_test = disp_corners - correction(corners)
    assert(np.allclose(correction_test[0], correction_test[1]))
    assert(np.allclose(correction_test[0], correction_test[2]))
    assert(np.allclose(correction_test[0], correction_test[3]))

    # Compute the displacements for pos
    disp = periodic(pos) - correction(pos)

    return disp

def dipole(self,
           sizemults: Tuple,
           boxtilt: bool = True,
           numreplicas: int = 5,
           shift: Optional[npt.ArrayLike] = None,
           shiftindex: Optional[int] = None,
           shiftscale: bool = False,
           center: Optional[npt.ArrayLike] = None,
           centerscale: bool = False,
           return_base_system: bool = False
           ) -> Union[System, Tuple[System, System]]:
    """
    Constructs a dislocation dipole configuration as described by Li, Wang,
    Chang, Cai, Bulatov, Ho, Yip, Phys Rev B 70(10) (2004) 104113 
    https://doi.org/10.1103/Physrevb.70.104113.  Two parallel dislocations with
    opposite Burgers vectors will be inserted into the system using a Volterra
    solution.  The resulting configuration is periodic in all three directions
    and constitutes a regularly-spaced 2D grid of parallel dislocations.  A
    shear strain is also applied to the system of 1/2 the Burgers vector to
    counteract the elastic strain of the dislocations and ensure stability.

    Parameters
    ----------
    sizemults : tuple
        The three size multipliers to use when generating the system.  Values
        should be positive integers if boxtilt is False.  When boxtilt is True,
        the multipliers are limited to values that result in full lattice
        vectors once the tilt is added.  Depending on the system, fractional
        values may be possible, or some integer values not allowed.
    boxtilt : bool, optional
        If True (default) then a tilt will be applied to the system such that
        the resulting periodic configuration will be consistent with a
        "quadripole" representation in which each dislocation will be
        surrounded by dislocations of the opposite sign in both the m- and n-
        directions.  This is achieved by adding half of the box vector
        most aligned with the m-axis to the box vector most aligned with the
        n-axis.  A value of False will not tilt the system, so only the
        sizemults will be applied to the rotated cell. The non-tilted system
        will have dislocations of the same sign aligned along the n-axis.
    numreplicas : int, optional
        Indicates how many image cells are used for computing the displacement
        field of the dipole. A rectangular grid of images is used
        meaning that there will be NxN total cells evaluated, of which NxN-1
        will be image cells.  Default value is 5.
    shift : array-like object, optional
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
    center : array-like object or None, optional
        Indicates where the dislocations are positioned in the configuration
        relative to the default locations.  For dipole configurations, the
        default locations are at relative positions of (1/4, 1/2) and 
        (3/4, 1/2) of the box dimensions of the final configuration that
        correspond to the dislocation solution's m- and n-axes.
    centerscale : bool, optional
        If False (default), a given center value will be taken as absolute
        Cartesian.  If True, a given center will be taken relative to the
        rotated cell's box vectors.
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
    # Extract box vector orientations
    cutindex = self.cutindex
    #lineindex = self.lineindex
    motionindex = self.motionindex

    # Multiply primitive uvws by sizemults to get dipole rotation uvws
    uvws_dipole = (np.asarray(sizemults) * self.uvws_prim.T).T
    
    # Apply boxtilt
    if boxtilt:
        uvws_dipole[cutindex] += uvws_dipole[motionindex] / 2

    # Check that the uvws are lattice vectors
    if not np.allclose(uvws_dipole, np.asarray(np.rint(uvws_dipole), dtype='int64')):
        raise ValueError(f'sizemults (and boxtilt) did not result in int lattice vectors: {uvws_dipole}')

    # Create base_system
    base_system = self.ucell_prim.rotate(uvws_dipole)

    # Handle shift parameters
    if shift is not None or shiftindex is not None:
        self.set_shift(shift, shiftindex, shiftscale)
    shift = self.shift
    
    # Handle center parameter
    if center is None:
        center = np.array([0,0,0])
    else:
        center = np.asarray(center)
    if centerscale:
        center = self.rcell.box.vector_crystal_to_cartesian(center)

    # Apply the rigid body shift to the atoms
    base_system.atoms.pos += shift

    # Move the box origin to the center of the system
    neworigin = - (base_system.box.vects[cutindex] + base_system.box.vects[motionindex]) / 2
    base_system.box_set(vects=base_system.box.vects, origin=neworigin)

    # Identify x1 and x2 positions
    m = self.dislsol.m
    length = base_system.box.vects.dot(m).dot(m)
    x2 = (length / 4) * m
    x1 = -x2

    # Shift atoms so that dislocation solution is centered at x1
    base_system.atoms.pos += x1
    base_system.wrap()

    # Compute displacement solution
    mvect = base_system.box.vects[motionindex]
    nvect = base_system.box.vects[cutindex]
    disp = self.dipole_displacement(self.dislsol, base_system.atoms.pos - center,
                                    x1, x2, mvect, nvect, N=numreplicas)

    # Create the dislocation system
    disl_system = deepcopy(base_system)
    disl_system.atoms.pos += disp

    # Apply the balancing strain (ONLY FOR SCREW RIGHT NOW!)
    newvects = disl_system.box.vects
    newvects[cutindex] -= 0.5 * self.dislsol.burgers
    disl_system.box_set(vects=newvects, origin=disl_system.box.origin, scale=True)

    self.set_systems(base_system, disl_system)

    if return_base_system:
        return base_system, disl_system
    else:
        return disl_system
