# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
import atommantest.core.Atoms
import atommantest.core.Box
import atommantest.core.System

#http://pymatgen.org
try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False

def load(structure):
    """
    Convert a pymatgen.Structure into an atomman.System.
    
    Parameters
    ----------
    structure : pymatgen.Structure
        A pymatgen representation of a structure.
    
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    elements : list
        A list of unique elemental symbols to pair with system's atypes.
    """
    
    assert has_pmg, 'pymatgen not imported'
    
    # Get box/lattice information
    box = atommantest.core.Box(vects = structure.lattice.matrix)
    
    # Get element information
    all_elements =  np.array([str(symbol) for symbol in structure.species])
    elements, atype = np.unique(all_elements, return_inverse = True)
    
    # Get atomic information
    prop = structure.site_properties
    prop['atype'] = atype + 1
    prop['pos'] = structure.cart_coords
    atoms = atommantest.core.Atoms(prop=prop)
    
    # Build system
    system = atommantest.core.System(atoms=atoms, box=box)
    
    return system, list(elements)
    