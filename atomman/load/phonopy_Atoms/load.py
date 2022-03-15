# coding: utf-8

# Standard Python imports
from typing import Optional

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import Atoms, Box, System

# https://atztogo.github.io/phonopy/
try:
    from phonopy.structure.atoms import PhonopyAtoms
    has_phonopy = True
except:
    has_phonopy = False

def load(phonopyatoms: PhonopyAtoms,
         symbols: Optional[tuple] = None,
         prop: Optional[dict] = None) -> System:
    """
    Convert a phonopy.structure.atoms.PhonopyAtoms into an atomman.System.
    
    Parameters
    ----------
    phonopyatoms : phonopy.structure.atoms.PhonopyAtoms
        A phonopy representation of a collection of atoms, based on ase.Atoms.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags.
    prop : dict, optional
        Dictionary containing any extra per-atom properties to include.
        
    Returns
    -------
    system : atomman.System
        A atomman representation of a system.
    """
    
    assert has_phonopy, 'phonopy not imported'
    
    if prop is None:
        prop = {}

    # Get box/cell information
    box = Box(vects = phonopyatoms.get_cell())
    pbc = (True, True, True)
    
    # Get element information
    all_elements = np.array(phonopyatoms.get_chemical_symbols())
    elements, atype = np.unique(all_elements, return_inverse = True)
    if symbols is None:
        symbols = elements

    # Get atomic information
    prop['atype'] = atype + 1
    prop['pos'] = phonopyatoms.get_positions()
    atoms = Atoms(prop=prop)
    
    # Build system
    system = System(atoms=atoms, box=box, pbc=pbc, symbols=symbols)
    
    return system