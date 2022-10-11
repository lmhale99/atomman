# coding: utf-8

# Standard Python imports
from typing import Optional

# http://www.numpy.org/
import numpy as np

# http://pymatgen.org
try:
    # import from newer pymatgen
    from pymatgen.core import Lattice, Structure

except ModuleNotFoundError:
    from pymatgen import Lattice, Structure

def dump(system,
         symbols: Optional[tuple] = None
         ) -> Structure:
    """
    Convert an atomman.System into a pymatgen.Structure.
    
    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    symbols : tuple, optional
        List of the element symbols that correspond to the atom types.  If not
        given, will use system.symbols if set, otherwise no element content
        will be included.
    
    Returns
    -------
    structure : pymatgen.Structure
        A pymatgen representation of a structure.
    """
    
    # Get box/lattice information
    lattice = Lattice(system.box.vects)
    
    # Get symbols information
    if symbols is None:
        symbols = system.symbols
    symbols = np.asarray(symbols)
    if None in symbols:
        raise ValueError('Symbols needed for all atypes')
    
    # Convert short symbols list to full species list
    atype = system.atoms.atype
    species = symbols[atype-1]
    
    # Get atomic information
    sites = system.atoms_prop(key='pos', scale=True)
    prop = {}
    for p in system.atoms_prop():
        if p != 'atype' and p != 'pos':
            prop[p] = system.atoms_prop(key=p)
    
    # Build structure
    structure = Structure(lattice, species, sites, site_properties=prop)
    
    return structure