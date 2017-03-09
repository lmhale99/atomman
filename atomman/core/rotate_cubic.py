import numpy as np
from copy import deepcopy
from math import ceil
from collections import OrderedDict

import atomman as am

def rotate_cubic(system, axes):
    """
    Rotate a cubic system according to the specified crystallographic axes.
    Returns an orthogonal system whose box lengths are:
        a * sqrt(i^2+j^2+k^2)
    where a is the original cubic box length, and ijk are the crystallographic
    indicies for the specific direction.  
    """
    
    
    axes = np.asarray(axes)

    #rotate system generically
    system = am.rotate(system, axes)

    #test for cubic
    a = system.box.a
    try:
        assert np.isclose(a, system.box.b), str(a) + ' ' + str(system.box.b)
        assert np.isclose(a, system.box.c), str(a) + ' ' + str(system.box.c)
        assert np.isclose(90.0, system.box.alpha), str(system.box.alpha)
        assert np.isclose(90.0, system.box.beta), str(system.box.beta)
        assert np.isclose(90.0, system.box.gamma), str(system.box.gamma)
    except:
        raise ValueError('Cubic system not given')
    
    #Test for integer axes values
    try:
        for ax_val in axes.flat:
            assert np.isclose(ax_val, int(ax_val))
    except:
        raise ValueError('axes values must be integers')    
    
    #Get magnitudes of the axes
    mag = np.linalg.norm(axes, axis=1)
    
    #compute number of atoms for the new system
    natoms = system.natoms*mag[0]*mag[1]*mag[2]

    if np.isclose(int(round(natoms)), natoms):
        natoms = int(round(natoms))
    else:
        raise ValueError('not an integer number of atoms associated with the axes.')

    #create a new box with scaled lattice parameters
    box = am.Box(a=a*mag[0], b=a*mag[1], c=a*mag[2], origin=system.box.origin)
    
    #supersize the rotated box
    m = int(ceil(max(mag)))
    system = am.supersize(system, (-m, m), (-m, m), (-m, m))

    #filter atoms for only those in the new box
    data = system.atoms.data

    #round x-positions near xlo, xhi and only keep data for atoms where xlo <= x < xhi
    data[np.where(np.isclose(data.T[1], box.xlo, atol=1e-8, rtol=0)), 1] = box.xlo 
    data = data[data.T[1] >= box.xlo]
    data[np.where(np.isclose(data.T[1], box.xhi, atol=1e-8, rtol=0)), 1] = box.xhi
    data = data[data.T[1] < box.xhi]
    
    #round y-positions near ylo, yhi and only keep data for atoms where ylo <= y < yhi
    data[np.where(np.isclose(data.T[2], box.ylo, atol=1e-8, rtol=0)), 2] = box.ylo 
    data = data[data.T[2] >= box.ylo]
    data[np.where(np.isclose(data.T[2], box.yhi, atol=1e-8, rtol=0)), 2] = box.yhi
    data = data[data.T[2] < box.yhi]
    
    #round z-positions near zlo, zhi and only keep data for atoms where zlo <= z < zhi
    data[np.where(np.isclose(data.T[3], box.zlo, atol=1e-8, rtol=0)), 3] = box.zlo 
    data = data[data.T[3] >= box.zlo]
    data[np.where(np.isclose(data.T[3], box.zhi, atol=1e-8, rtol=0)), 3] = box.zhi
    data = data[data.T[3] < box.zhi]
    
    #deepcopy the data array to guarantee that it is new and separate
    data = deepcopy(data)
    

    #rebuild views
    start = 0
    view = OrderedDict()
    for k in system.atoms.view:
        vshape = (len(data), ) + system.atoms.view[k].shape[1:]
        view[k] = data[:, start : start + system.atoms.view[k][0].size]
        view[k].shape = vshape
        start = start + system.atoms.view[k][0].size

    #create atoms from natoms, data, view and dtype
    atoms = am.Atoms(data=data, view=view, prop_dtype=system.atoms.dtype)
    
    assert natoms == atoms.natoms, 'natom mismatch after rotation!'
    
    #return new system
    return am.System(box=box, atoms=atoms)