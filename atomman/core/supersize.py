from .Atoms import Atoms
from .Box import Box
from .System import System

import numpy as np


def supersize(system, a_size, b_size, c_size):
    """
    Builds a large system based on a seed system and multipliers.
    
    Keyword Arguments:
    system -- atomman.System to use as the seed.
    a_size -- int or tuple of 2 ints for multiplying along the avect direction.
    b_size -- int or tuple of 2 ints for multiplying along the bvect direction.
    c_size -- int or tuple of 2 ints for multiplying along the cvect direction.
    
    The multiplier values *_size are taken to be integer tuples (m, n) where m <= 0 and n >= 0.
    The system multiplication works such that if n = -m, then the seed system's origin will be at the center of the new system.
    If only one integer is given, then it is assigned to m or n depending on its sign, and the other value is taken to be 0.  
    """
    #initial parameter setup
    sizes = [a_size, b_size, c_size]
    mults = np.array([0, 0, 0], dtype=int)
    vects = system.box.vects
    origin = system.box.origin
    spos = system.atoms_prop(key='pos', scale=True)
    
    
    for i in xrange(3):
        #check values in sizes
        if isinstance(sizes[i], (int, long)):
            if sizes[i] > 0:
                sizes[i] = (0, sizes[i])
            elif sizes[i] < 0:
                sizes[i] = (sizes[i], 0)
        elif isinstance(sizes[i], tuple):
            assert len(sizes[i]) == 2, 'Invalid system multipliers'
            assert isinstance(sizes[i][0], (int, long)) and sizes[i][0] <= 0, 'Invalid system multipliers'
            assert isinstance(sizes[i][1], (int, long)) and sizes[i][1] >= 0, 'Invalid system multipliers'
        else:
            raise TypeError('Invalid system multipliers')
        
        #calculate multipliers and scale box and first set of positions accordingly
        mults[i] = sizes[i][1] - sizes[i][0]
        assert mults[i] != 0, 'Cannot multiply system dimension by zero'
        spos[:,i] /= mults[i]
        origin += vects[i] * sizes[i][0]
        vects[i] *= mults[i]
        
    #initilize new Box and Atoms
    box = Box(vects=vects, origin=origin)
    natoms = system.natoms * mults[0] * mults[1] * mults[2]
    atoms = Atoms(natoms=natoms)
    
    #Copy over all property values using numpy broadcasting
    for prop in system.atoms_prop():
        o_values = system.atoms_prop(key=prop)
        n_values = np.empty((mults[0] * mults[1] * mults[2],) + o_values.shape, dtype = o_values.dtype)
        n_values[:] = system.atoms_prop(key=prop)
        n_shape = n_values.shape
        n_shape = (n_shape[0]*n_shape[1],) + n_shape[2:]
        atoms.prop(key=prop, value=n_values.reshape(n_shape))
    
    #Expand spos using broadcasting
    n_spos = np.empty((mults[0] * mults[1] * mults[2],) + spos.shape)
    n_spos[:] = spos
    n_shape = n_spos.shape
    n_shape = (n_shape[0]*n_shape[1],) + n_shape[2:]
    n_spos = n_spos.reshape(n_shape)
    
    #use broadcasting to create arrays to add to spos
    test = np.empty(mults[0]*system.natoms)
    test.shape = (system.natoms, mults[0])
    test[:] = np.arange(mults[0])
    x = test.T.flatten()

    test = np.empty(mults[1]*len(x))
    test.shape = (len(x), mults[1])
    test[:] = np.arange(mults[1])
    y = test.T.flatten()
    test.shape = (mults[1], len(x))
    test[:] = x
    x = test.flatten()

    test = np.empty(mults[2]*len(x))
    test.shape = (len(x), mults[2])
    test[:] = np.arange(mults[2])
    z = test.T.flatten()
    test.shape = (mults[2], len(x))
    test[:] = x
    x = test.flatten()
    test[:] = y
    y = test.flatten()
    
    #xyz is displacement values to add to spos
    xyz = np.hstack((x[:,np.newaxis], y[:,np.newaxis], z[:,np.newaxis])) * np.array([1./mults[0], 1./mults[1], 1./mults[2]])
    
    #save pos values, return new System
    atoms.prop(key='pos', value=n_spos + xyz)
    
    return System(box=box, atoms=atoms, scale=True)
        
        
    
        
    
    
    
    