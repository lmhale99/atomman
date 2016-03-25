import numpy as np
import atomman as am
from copy import deepcopy
from math import ceil
from collections import OrderedDict

def rotate(system, axes):
    """
    Rotate a system by transforming both atomic positions and lattice vectors
    according to axes. 
    
    Arguments:
    system -- atomman.System that is being rotated.
    axes -- three right-handed orthogonal vectors defining the transformation.
            The standard transformation matrix is taken by normalizing axes.
    
    The system returned should be identical to the original
    system, just with a new choice of axes.    
    """
    #normalize axes to get the transformation matrix, T
    T = am.tools.axes_check(axes)
    
    #transform the box vectors using T
    avect = T.dot(system.box.avect)
    bvect = T.dot(system.box.bvect)
    cvect = T.dot(system.box.cvect)
    
    #create a new system using deepcopy
    new_system = deepcopy(system)
    
    #transform the system's axes
    new_system.box_set(avect=avect, bvect=bvect, cvect=cvect, scale=True)
    
    return new_system
    

def rotate_cubic(system, axes):
    """
    Rotate a cubic system according to the specified crystallographic axes.
    Returns an orthogonal system whose box lengths are:
        a * sqrt(i^2+j^2+k^2)
    where a is the original cubic box length, and ijk are the crystallographic
    indicies for the specific direction.  
    """
    #rotate system generically
    system = rotate(system, axes)
    
    #test for cubic
    a = system.box.a
    if True:
    #try:
        assert np.isclose(a, system.box.b), str(a) + ' ' + str(system.box.b)
        assert np.isclose(a, system.box.c), str(a) + ' ' + str(system.box.c)
        assert np.isclose(90.0, system.box.alpha), str(system.box.alpha)
        assert np.isclose(90.0, system.box.beta), str(system.box.beta)
        assert np.isclose(90.0, system.box.gamma), str(system.box.gamma)
    else:
    #except:
        raise ValueError('Cubic system not given')
    
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
    system = am.tools.supersize(system, (-m, m), (-m, m), (-m, m))

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
    atoms = am.Atoms(natoms=natoms, data=data, view=view, prop_dtype=system.atoms.dtype)

    #return new system
    return am.System(box=box, atoms=atoms)