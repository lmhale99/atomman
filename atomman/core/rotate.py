import atomman as am
from copy import deepcopy

from ..tools import axes_check

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
    T = axes_check(axes)
    
    #transform the box vectors using T
    avect = T.dot(system.box.avect)
    bvect = T.dot(system.box.bvect)
    cvect = T.dot(system.box.cvect)
    
    #create a new system using deepcopy
    new_system = deepcopy(system)
    
    #transform the system's axes
    new_system.box_set(avect=avect, bvect=bvect, cvect=cvect, scale=True)
    
    return new_system
    

