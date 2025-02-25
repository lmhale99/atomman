# coding: utf-8
# Standard Python libraries
from typing import Optional

# https://atztogo.github.io/spglib/python-spglib.html
import spglib

# atomman imports
from ... import System, load

def dump(system: System,
         symprec: float = 1e-5,
         normalize: Optional[str] = 'lammps') -> System:
    """
    Wrapper around spglib.find_primitive() that finds a primitive unit
    cell for the system.  This works best for small systems where atoms are
    already close to or at ideal positions.

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
    
    # Convert to spglib cell
    cell = system.dump('spglib_cell')
    
    # Use spglib to identify the primitive cell
    primitive_cell = spglib.find_primitive(cell, symprec=symprec)
    
    # Return original system if primitive search failed
    if primitive_cell is None:
        return system
    
    # Convert back to an atomman System and normalize
    primitive_system = load('spglib_cell', primitive_cell, symbols=system.symbols)
    
    # Normalize
    if normalize is not None:
        primitive_system = primitive_system.normalize(normalize)
        
    return primitive_system