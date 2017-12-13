# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
import atomman.core.Atoms
import atomman.core.Box
import atomman.core.System

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
    box = atomman.core.Box(vects = structure.lattice.matrix)
    
    # Get element information
    all_elements =  np.array([str(symbol) for symbol in structure.species])
    elements, atype = np.unique(all_elements, return_inverse = True)
    
    # Get atomic information
    prop = structure.site_properties
    prop['atype'] = atype + 1
    prop['pos'] = structure.cart_coords
    atoms = atomman.core.Atoms(prop=prop)
    
    # Build system
    system = atomman.core.System(atoms=atoms, box=box)
    
    return system, list(elements)
    
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