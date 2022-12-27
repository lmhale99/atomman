# coding: utf-8
from typing import Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
from . import VolterraDislocation

def dislocation_dipole_displacement(dislsol: VolterraDislocation,
                                    pos: npt.ArrayLike,
                                    x1: Union[float, npt.ArrayLike],
                                    x2: Union[float, npt.ArrayLike],
                                    mvect: npt.ArrayLike,
                                    nvect: npt.ArrayLike,
                                    N: int = 5) -> np.ndarray:
    """
    Uses the method of Cai, Bulatov, Chang, Li & Yip, Phil Mag 2003, 83(5), 539–567
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
    N : int
        Indicates how many image cells are used. A rectangular grid of images is used
        meaning that there will be NxN total cells evaluated, of which NxN-1 will be
        image cells.
    
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
    delta_disp
    
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