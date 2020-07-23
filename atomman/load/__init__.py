# coding: utf-8

class FileFormatError(Exception):
    pass

from .ase_Atoms import load as load_ase_Atoms
from .cif import load as load_cif
from .pymatgen_Structure import load as load_pymatgen_Structure
from .table import load as load_table
from .system_model import load as load_system_model
from .prototype import load as load_prototype
from .crystal import load as load_crystal
from .poscar import load as load_poscar
from .atom_data import load as load_atom_data
from .atom_dump import load as load_atom_dump
from .spglib_cell import load as load_spglib_cell
from .phonopy_Atoms import load as load_phonopy_Atoms

__all__ = ['FileFormatError', 'load', 'load_ase_Atoms', 'load_pymatgen_Structure',
           'load_table', 'load_system_model', 'load_poscar', 'load_atom_data',
           'load_atom_dump', 'load_cif', 'load_spglib_cell', 'load_phonopy_Atoms',
           'load_prototype', 'load_crystal']

def load(style, *args, **kwargs):
    """
    Load a System.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to load as an atomman.System
    args 
        Any positional-dependent arguments to pass to the underlying load methods.
    kwargs
        Any keyword arguments to pass to the underlying load methods.
        
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    """
    
    if style == 'system_model':
        return load_system_model(*args, **kwargs)
    
    elif style == 'prototype':
        return load_prototype(*args, **kwargs)

    elif style == 'crystal':
        return load_crystal(*args, **kwargs)

    elif style == 'cif':
        return load_cif(*args, **kwargs)
    
    elif style == 'atom_data':
        return load_atom_data(*args, **kwargs)
    
    elif style == 'atom_dump':
        return load_atom_dump(*args, **kwargs)
    
    elif style == 'table':
        return load_table(*args, **kwargs)
    
    elif style == 'ase_Atoms':
        return load_ase_Atoms(*args, **kwargs)
    
    elif style == 'phonopy_Atoms':
        return load_phonopy_Atoms(*args, **kwargs)
    
    elif style == 'pymatgen_Structure':
        return load_pymatgen_Structure(*args, **kwargs)
    
    elif style == 'poscar':
        return load_poscar(*args, **kwargs)
    
    elif style == 'spglib_cell':
        return load_spglib_cell(*args, **kwargs)
    
    else:
        raise ValueError('Unsupported style')