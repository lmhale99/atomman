import numpy as np
from dvect import dvect

def displacement(system_0, system_1):
    """Compute the displacement vectors between all matching atoms for two systems."""
    
    assert system_0.natoms == system_1.natoms,  'systems have different number of atoms'
    
    disp = dvect(system_0.atoms.view['pos'], system_1.atoms.view['pos'], system_1.box, system_1.pbc)
                            
    return disp
    