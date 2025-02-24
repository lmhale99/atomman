# coding: utf-8

# Standard Python imports
from typing import Tuple, Union

# http://www.numpy.org/
import numpy as np

# https://wiki.fysik.dtu.dk/ase/
import ase

def dump(system,
         symbols: Union[str, list, None] = None,
         return_prop: bool = False
         ) -> Union[ase.Atoms, Tuple[ase.Atoms, dict]]:
    """
    Convert an atomman.System into an ase.Atoms.
    
    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    symbols : tuple, optional
        List of the element symbols that correspond to the atom types.  If not
        given, will use system.symbols if set, otherwise no element content
        will be included.
    return_prop : bool, optional
        Indicates if the extra per-atom properties are to be returned in a
        dictionary.  Default value is False.
    Returns
    -------
    aseatoms : ase.Atoms
        An ase representation of a collection of atoms.
    prop : dict
        Dictionary containing any extra per-atom properties to include.
        Returned if return_prop is True.
    """
    
    # Get box/cell information
    cell = system.box.vects
    pbc = system.pbc
    
    # Get symbols information
    if symbols is None:
        symbols = system.symbols
    symbols = np.asarray(symbols)
    if None in symbols:
        raise ValueError('Symbols needed for all atypes')
    
    # Convert short symbols list to full allsymbols list
    atype = system.atoms.atype
    allsymbols = symbols[atype-1]
    
    # Get atomic information
    spos = system.atoms_prop('pos', scale=True)
    prop = {}
    for p in system.atoms_prop():
        if p != 'atype' and p != 'pos':
            prop[p] = system.atoms_prop(key=p)
    
    # Build Atoms
    aseatoms = ase.Atoms(symbols=allsymbols, scaled_positions=spos,
                         pbc=pbc, cell=cell)
    
    if return_prop is True:
        return aseatoms, prop
    else:
        return aseatoms