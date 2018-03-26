# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

#http://pymatgen.org
try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False

def dump(system, elements):
    """
    Convert an atomman.System into a pymatgen.Structure.
    
    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    elements : list
        A list of unique elemental symbols to pair with system's atypes.
    
    Returns
    -------
    structure : pymatgen.Structure
        A pymatgen representation of a structure.
    """
    
    assert has_pmg, 'pymatgen not imported'
    
    # Get box/lattice information
    lattice = pmg.Lattice(system.box.vects)
    
    # Get element/species information
    elements = np.asarray(elements)
    atype = system.atoms.atype
    species = elements[atype-1]
    
    # Get atomic information
    sites = system.atoms_prop(key='pos', scale=True)
    prop = {}
    for p in system.atoms_prop():
        if p != 'atype' and p != 'pos':
            prop[p] = system.atoms_prop(key=p)
    
    # Build structure
    structure = pmg.Structure(lattice, species, sites, site_properties=prop)
    
    return structure