# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
import atomman.core.Atoms
import atomman.core.Box
import atomman.core.System

# https://wiki.fysik.dtu.dk/ase/
try:
    import ase
    has_ase = True
except:
    has_ase = False

def load(aseatoms, prop={}):
    """
    Convert an ase.Atoms into an atomman.System and list of elements.
    
    Parameters
    ----------
    aseatoms : ase.Atoms
        An ase representation of a collection of atoms.
    prop : dict, optional
        Dictionary containing any extra per-atom properties to include.
        
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    elements : list
        A list of unique elemental symbols to pair with system's atypes.
    """
    
    assert has_ase, 'ase not imported'
    
    # Get box/cell information
    box = atomman.core.Box(vects = aseatoms.get_cell())
    pbc = aseatoms.get_pbc()
    
    # Get element information
    all_elements = np.array(aseatoms.get_chemical_symbols())
    elements, atype = np.unique(all_elements, return_inverse = True)
    
    # Get atomic information
    prop['atype'] = atype + 1
    prop['pos'] = aseatoms.get_positions()
    atoms = atomman.core.Atoms(prop=prop)
    
    # Build system
    system = atomman.core.System(atoms=atoms, box=box, pbc=pbc)
    
    return system, elements