import numpy as np

import atomman as am

def disregistry(basesystem, dislsystem):
    """
    Computes the disregistry profile for a dislocation system.  Assumes that
    the dislocation line is along the z-axis and the slip plane is the y=0
    plane.  These assumptions match how atomman.defect.Stroh generates
    dislocations.
    
    Parameters
    ----------
    basesystem : atomman.System
        A perfect reference system with atoms directly corresponding to atoms
        in dislsystem.
    dislsystem : atomman.System
        A dislocation-containing system.
        
    Returns
    -------
    x : numpy.ndarray
        The (N,) array of unique x-coordinates (atomic columns) neighboring
        the slip plane.
    disregistry : numpy.ndarray
        A (N, 3) array of the dislocation's disregistry at each x.
    """
    # Extract pos from basesystem and compute atomic displacements
    basepos = basesystem.atoms.view['pos'] 
    displacement = am.displacement(basesystem, dislsystem)
    
    # Identify atoms just above and below slip plane
    uniquey = np.unique(basepos[:, 1])
    abovey = uniquey[uniquey > 0].min()
    belowy = uniquey[uniquey < 0].max()
    
    # Get x-coordinates of atoms just above and below slip plane
    abovex = basepos[np.isclose(basepos[:, 1], abovey)][:, 0]
    belowx = basepos[np.isclose(basepos[:, 1], belowy)][:, 0]
    
    # Identify unique x-coordinates
    uabovex = np.unique(abovex)
    ubelowx = np.unique(belowx)
    x = np.union1d(uabovex, ubelowx)
    
    # Get displacements of atoms just above and below slip plane
    abovedisp = displacement[np.isclose(basepos[:, 1], abovey)]
    belowdisp = displacement[np.isclose(basepos[:, 1], belowy)]
    
    # Average displacements for the same x-coordinates
    abovedispmean = np.empty((len(uabovex), 3))
    for i, ix in enumerate(uabovex):
        abovedispmean[i] = abovedisp[np.isclose(abovex, ix)].mean(axis=0)
    belowdispmean = np.empty((len(ubelowx), 3))
    for i, ix in enumerate(ubelowx):
        belowdispmean[i] = belowdisp[np.isclose(belowx, ix)].mean(axis=0)
    
    # Linearly interpolate displacement values to all x
    abovedispinterp = np.vstack([np.interp(x, uabovex, abovedispmean[:,0]),
                                 np.interp(x, uabovex, abovedispmean[:,1]),
                                 np.interp(x, uabovex, abovedispmean[:,2])]).T
    belowdispinterp = np.vstack([np.interp(x, ubelowx, belowdispmean[:,0]),
                                 np.interp(x, ubelowx, belowdispmean[:,1]),
                                 np.interp(x, ubelowx, belowdispmean[:,2])]).T
    
    # Calculate disregistry
    disregistry = abovedispinterp - belowdispinterp
    
    return x, disregistry