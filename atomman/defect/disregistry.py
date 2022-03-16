# coding: utf-8

# Standard Python libraries
from typing import Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from .. import displacement, System

def disregistry(basesystem: System,
                dislsystem: System,
                m: npt.ArrayLike = [1.0, 0.0, 0.0],
                n: npt.ArrayLike = [0.0, 1.0, 0.0],
                planepos: npt.ArrayLike = [0.0, 0.0, 0.0]
                ) ->  Tuple[np.ndarray, np.ndarray]:
    """
    Computes the disregistry profile for a dislocation system.
    
    Parameters
    ----------
    basesystem : atomman.System
        A perfect reference system with atoms directly corresponding to atoms
        in dislsystem.
    dislsystem : atomman.System
        A dislocation-containing system.
    m : array-like object, optional
        The dislocation solution m unit vector.  This vector is in the slip
        plane and perpendicular to the dislocation line direction.  Default
        value is [1,0,0] (Cartesian x-axis).
    n : array-like object, optional
        The dislocation solution n unit vector.  This vector is normal to the 
        slip plane.  Only needed if dislsol is not given.  Default
        value is [0,1,0] (Cartesian y-axis).
    planepos : array-like object, optional
        A position on the slip plane so that the plane can be fully defined.
        The slip plane position should fall between two planes of atoms.
        Default value is [0,0,0].
    
    Returns
    -------
    coord : numpy.ndarray
        The (N,) array of unique coord-coordinates (atomic columns) neighboring
        the slip plane.
    disregistry : numpy.ndarray
        A (N, 3) array of the dislocation's disregistry at each coord.
    """
    # Handle m, n and planepos
    m = np.asarray(m, dtype=float)
    n = np.asarray(n, dtype=float)
    planepos = np.asarray(planepos, dtype=float)
    #if not np.isclose(m.dot(planepos), 0.0, atol=1e-8, rtol=0.0):
    #    raise ValueError('m and planepos must be perpendicular')
    
    # Extract pos from basesystem and compute atomic displacements
    basepos = basesystem.atoms.pos
    disp = displacement(basesystem, dislsystem)
    
    # Transform basepos and planepos to calculation coordinates
    allx = np.dot(basepos, m)
    ally = np.dot(basepos, n)
    midy = np.dot(planepos, n)
    
    # Identify atoms just above and below slip plane
    uniquey = np.unique(ally)
    abovey = uniquey[uniquey > midy].min()
    belowy = uniquey[uniquey < midy].max()
    if np.isclose(abovey, belowy):
        raise ValueError('planepos must fall between atomic planes')
    #uniquey = np.unique(basepos[:, 1])
    #abovey = uniquey[uniquey > 0].min()
    #belowy = uniquey[uniquey < 0].max()
    
    # Get coordinates of atoms just above and below slip plane
    abovex = allx[np.isclose(ally, abovey)]
    belowx = allx[np.isclose(ally, belowy)]
    #abovex = basepos[np.isclose(basepos[:, 1], abovey)][:, 0]
    #belowx = basepos[np.isclose(basepos[:, 1], belowy)][:, 0]
    
    # Identify unique coordinates
    uabovex = np.unique(abovex)
    ubelowx = np.unique(belowx)
    coord = np.union1d(uabovex, ubelowx)
    
    # Get displacements of atoms just above and below slip plane
    abovedisp = disp[np.isclose(ally, abovey)]
    belowdisp = disp[np.isclose(ally, belowy)]
    #abovedisp = disp[np.isclose(basepos[:, 1], abovey)]
    #belowdisp = disp[np.isclose(basepos[:, 1], belowy)]
    
    # Average displacements for the same coordinates
    abovedispmean = np.empty((len(uabovex), 3))
    for i, ix in enumerate(uabovex):
        abovedispmean[i] = abovedisp[np.isclose(abovex, ix)].mean(axis=0)
    belowdispmean = np.empty((len(ubelowx), 3))
    for i, ix in enumerate(ubelowx):
        belowdispmean[i] = belowdisp[np.isclose(belowx, ix)].mean(axis=0)
    
    # Linearly interpolate displacement values to all coord
    abovedispinterp = np.vstack([np.interp(coord, uabovex, abovedispmean[:,0]),
                                 np.interp(coord, uabovex, abovedispmean[:,1]),
                                 np.interp(coord, uabovex, abovedispmean[:,2])]).T
    belowdispinterp = np.vstack([np.interp(coord, ubelowx, belowdispmean[:,0]),
                                 np.interp(coord, ubelowx, belowdispmean[:,1]),
                                 np.interp(coord, ubelowx, belowdispmean[:,2])]).T
    
    # Calculate disregistry
    disregistry = abovedispinterp - belowdispinterp
    
    return coord, disregistry