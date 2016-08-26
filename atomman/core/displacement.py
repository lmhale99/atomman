import numpy as np
from dvect import dvect

def displacement(system_0, system_1, box_reference='final'):
    """Compute the displacement vectors between all matching atoms for two systems."""
    
    assert system_0.natoms == system_1.natoms,  'systems have different number of atoms'
    
    if box_reference == 'final':
        disp = dvect(system_0.atoms.view['pos'], system_1.atoms.view['pos'], system_1.box, system_1.pbc)
    elif box_reference == 'initial':
        disp = dvect(system_0.atoms.view['pos'], system_1.atoms.view['pos'], system_0.box, system_0.pbc)
    elif box_reference is None:
        disp = system_1.atoms.view['pos'] - system_0.atoms.view['pos']
    else:
        raise ValueError("box_reference must be 'final', 'initial', or None")
    
    return disp
    