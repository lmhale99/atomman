# coding: utf-8
# Standard Python libraries
from typing import Optional

# https://atztogo.github.io/spglib/python-spglib.html
try:
    import spglib
    has_spglib = True
except ImportError:
    has_spglib = False


# atomman imports
from ... import System, load_spglib_cell

def dump(system: System,
         symprec: float = 1e-5,
         normalize: Optional[str] = 'lammps') -> System:
    """
    Converts a given system into a primitive unit cell.  NOTE: This works best for small systems,
    e.g. conventional unit cells or relatively small and standard supercells.

    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    symprec : float, optional
        Absolute length tolerance to use in identifying symmetry of atomic
        sites and system boundaries. Default value is 1e-5.
    normalize : str or None
        Indicates which normalization scheme, if any, to use on the identified
        primitive cell.  None will return exactly as obtained from spglib.
        Default value is 'lammps', meaning that the cell will be compatible
        with LAMMPS.
    """
    assert has_spglib, 'spglib must be installed to use dump primitive cell'

    # Convert to spglib cell
    cell = system.dump('spglib_cell')
    
    # Use spglib to identify the primitive cell
    primitive_cell = spglib.find_primitive(cell, symprec=symprec)
    
    # Return original system if primitive search failed
    if primitive_cell is None:
        return system
    
    # Convert back to an atomman System and normalize
    primitive_system = load_spglib_cell(primitive_cell, symbols=system.symbols)
    
    # Normalize
    if normalize is not None:
        primitive_system = primitive_system.normalize(normalize)
        
    return primitive_system