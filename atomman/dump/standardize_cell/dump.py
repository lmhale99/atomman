# coding: utf-8
# Standard Python libraries
from typing import Optional

# https://atztogo.github.io/spglib/python-spglib.html
import spglib

# atomman imports
from ... import System, load

def dump(system: System,
         symprec: float = 1e-5,
         to_primitive: bool = False,
         no_idealize: bool = False,
         angle_tolerance: float = -1.0,
         normalize: Optional[str] = 'lammps') -> System:
    """
    Wrapper around spglib.standardize_cell() that finds a standardized unit
    cell for the system.  This works best for small systems where atoms are
    already close to or at ideal positions.

    Parameters
    ----------
    system : atomman.System
        An atomman representation of a system.
    symprec : float, optional
        Absolute length tolerance to use in identifying symmetry of atomic
        sites and system boundaries. Default value is 1e-5.
    to_primitive : bool, optional
        If True, the returned unit cell will be primitive.  If False (default),
        the returned unit cell will be conventional.
    no_idealize : bool, optional
        If False (default) the lengths and angles of box vectors and atomic
        positions will all be idealized based on symmetry.  Setting this to
        True turns off the idealization.
    angle_tolerance : float, optional
        Symmetry search tolerance in the unit of angle deg.  The spglib
        documentation suggests that this is normally not to be used, and giving
        a negative value (default) will use an optimized internal algorithm
        instead.
    normalize : str or None
        Indicates which normalization scheme, if any, to use on the identified
        primitive cell.  None will return exactly as obtained from spglib.
        Default value is 'lammps', meaning that the cell will be compatible
        with LAMMPS.
    """
    
    # Convert to spglib cell
    cell = system.dump('spglib_cell')
    
    # Use spglib to identify the primitive cell
    standardized_cell = spglib.standardize_cell(cell, symprec=symprec,
                                                to_primitive=to_primitive,
                                                no_idealize=no_idealize,
                                                angle_tolerance=angle_tolerance)
    
    # Return original system if standardization failed
    if standardized_cell is None:
        return system
    
    # Convert back to an atomman System and normalize
    standardized_system = load('spglib_cell', standardized_cell, symbols=system.symbols)
    
    # Normalize
    if normalize is not None:
        standardized_system = standardized_system.normalize(normalize)
        
    return standardized_system