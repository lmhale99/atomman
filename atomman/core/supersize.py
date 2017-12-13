# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from . import Atoms
from . import Box
import atomman.core.System
from ..compatibility import range, inttype

def supersize(system, a_size, b_size, c_size):
    """
    Creates a larger system from a given system by replicating it along the
    system's box vectors.
    
    The multiplier values *_size are taken to be integer tuples (m, n) where
    m <= 0 and n >= 0.  The system multiplication works such that if n = -m,
    then the seed system's origin will be at the center of the new system.
    If only one integer is given, then it is assigned to m or n depending on
    its sign, and the other value is taken to be 0.
    
    Parameters
    ----------
    system : atomman.System 
        The seed system to use in generating the larger system.
    a_size : int or tuple of int
        Single int or two integers specifying replication along the avect
        direction.
    b_size -- int or tuple of int
        Single int or two integers specifying replication along the bvect
        direction.
    c_size -- int or tuple of int
        Single int or two integers specifying replication along the cvect
        direction.
    
    Returns
    -------
    atomman.System
        A new system created by replicating the given seed system according to
        the *_size parameters.
      
    """
    # Extract parameters
    sizes = [a_size, b_size, c_size]
    mults = np.array([0, 0, 0], dtype=int)
    vects = system.box.vects
    origin = system.box.origin
    spos = system.atoms_prop('pos', scale=True)
    
    # Check the *_size values
    for i in range(3):
        
        # Change single int to tuple of two int
        if isinstance(sizes[i], inttype):
            if sizes[i] > 0:
                sizes[i] = (0, sizes[i])
            elif sizes[i] < 0:
                sizes[i] = (sizes[i], 0)
        
        elif isinstance(sizes[i], tuple):
            if True:
                assert len(sizes[i]) == 2, str(len(sizes[i]))
                assert isinstance(sizes[i][0], inttype), str(sizes[i][0])
                assert sizes[i][0] <= 0, str(sizes[i][0])
                assert isinstance(sizes[i][1], inttype), str(sizes[i][1])
                assert sizes[i][1] >= 0, str(sizes[i][1])
            else:
                raise TypeError('Invalid system multipliers')
        else:
            raise TypeError('Invalid system multipliers')
        
        # Calculate full multipliers 
        mults[i] = sizes[i][1] - sizes[i][0]
        if mults[i] == 0:
            raise ValueError('Cannot multiply system dimension by zero')
            
        # Scale box and first set of positions accordingly
        spos[:,i] /= mults[i]
        origin += vects[i] * sizes[i][0]
        vects[i] *= mults[i]
        
    # Initialize new Box and Atoms
    box = Box(vects=vects, origin=origin)
    natoms = system.natoms * mults[0] * mults[1] * mults[2]
    atoms = Atoms(natoms=natoms)
    
    # Copy over all property values (except pos) using numpy broadcasting
    for key in system.atoms_prop():
        if key == 'pos':
            continue
    
        # Get old array
        old = system.atoms.view[key]
        
        # Create new array and broadcast old to it
        new = np.empty((mults[0] * mults[1] * mults[2],) + old.shape, dtype = old.dtype)
        new[:] = old
        
        # Reshape new and save to atoms
        new_shape = new.shape
        new_shape = (new_shape[0] * new_shape[1],) + new_shape[2:]
        atoms.view[key] = np.array(new.reshape(new_shape))
    
    # Expand spos using broadcasting
    new_spos = np.empty((mults[0] * mults[1] * mults[2],) + spos.shape)
    new_spos[:] = spos
    
    # Reshape spos
    new_shape = new_spos.shape
    new_shape = (new_shape[0]*new_shape[1],) + new_shape[2:]
    new_spos = new_spos.reshape(new_shape)
    
    # Use broadcasting to create arrays to add to spos
    test = np.empty(mults[0] * system.natoms)
    test.shape = (system.natoms, mults[0])
    test[:] = np.arange(mults[0])
    x = test.T.flatten()

    test = np.empty(mults[1] * len(x))
    test.shape = (len(x), mults[1])
    test[:] = np.arange(mults[1])
    y = test.T.flatten()
    test.shape = (mults[1], len(x))
    test[:] = x
    x = test.flatten()

    test = np.empty(mults[2] * len(x))
    test.shape = (len(x), mults[2])
    test[:] = np.arange(mults[2])
    z = test.T.flatten()
    test.shape = (mults[2], len(x))
    test[:] = x
    x = test.flatten()
    test[:] = y
    y = test.flatten()
    
    # xyz is displacement values to add to spos
    xyz = (np.hstack((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis])) 
         * np.array([1 / mults[0], 1 / mults[1], 1 / mults[2]]))
    
    # Save pos values, return new System
    atoms.view['pos'] = new_spos + xyz
    
    return atomman.core.System(box=box, atoms=atoms, scale=True)