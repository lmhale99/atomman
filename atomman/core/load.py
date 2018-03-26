# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

import atomman.convert
from ..compatibility import range

def load(style, input, **kwargs):
    """
    Load a System.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to load as an atomman.System
    input : str, file-like object or object
        The content to load.
    kwargs
        Any extra keyword arguments to pass to the underlying load methods.
        
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    symbols : list
        The list of atomic symbols corresponding to the system's atom types.
        Will be a list of None if symbol information is not in input.
    """
    
    if style == 'system_model':
        system, symbols = atomman.convert.system_model.load(input, **kwargs)
    
    elif style == 'cif':
        system, symbols = atomman.convert.cif.load(input, **kwargs)
    
    elif style == 'atom_data':
        system = atomman.convert.atom_data.load(input, **kwargs)
        symbols = [None for i in range(system.natypes)]
    
    elif style == 'atom_dump':
        system = atomman.convert.atom_dump.load(input, **kwargs)
        symbols = [None for i in range(system.natypes)]
    
    elif style == 'table':
        system = atomman.convert.table.load(input, **kwargs)
        symbols = [None for i in range(system.natypes)]
    
    elif style == 'ase_Atoms':
        system, symbols = atomman.convert.ase_Atoms.load(input, **kwargs)
    
    elif style == 'pymatgen_Structure':
        system, symbols = atomman.convert.pymatgen_Structure.load(input, **kwargs)
    
    elif style == 'poscar':
        system, symbols = atomman.convert.poscar.load(input, **kwargs)
    
    else:
        raise ValueError('Unsupported style')
    
    return system, symbols