# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

import atommantest.convert

def dump(style, system, f=None, **kwargs):
    """
    Convert a System to another format.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to dump the atomman.System as.
    system : atomman.System
        The system to convert.
    f : str or file-like object
        File path or buffer to save content to for the text-based formats.
    kwargs
        Any extra keyword arguments to pass to the underlying dump methods.
        
    Returns
    -------
    str, object or tuple
        Any content returned by the underlying dump methods.
    
    """
    
    if style == 'system_model':
        return atommantest.convert.system_model.dump(system, f=f, **kwargs)
    
    elif style == 'atom_data':
        return atommantest.convert.atom_data.dump(system, f=f, **kwargs)
    
    elif style == 'atom_dump':
        return atommantest.convert.atom_dump.dump(system, f=f, **kwargs)
    
    elif style == 'table':
        return atommantest.convert.table.dump(system, f=f, **kwargs)
    
    elif style == 'ase_Atoms':
        assert f is None, 'f not allowed with style ase_Atoms'
        return atommantest.convert.ase_Atoms.dump(system, **kwargs)
    
    elif style == 'pymatgen_Structure':
        assert f is None, 'f not allowed with style pymatgen_Structure'
        return atommantest.convert.pymatgen_Structure.dump(system, **kwargs)
    
    elif style == 'poscar':
        return atommantest.convert.poscar.dump(system, f=f, **kwargs)
    
    else:
        raise ValueError('Unsupported style')