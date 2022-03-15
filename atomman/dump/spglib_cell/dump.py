# coding: utf-8

def dump(system) -> tuple:
    """
    Convert an atomman.System into the cell tuple used by spglib.
    
    Parameters
    ----------
    system : atomman.System
        A atomman representation of a system.
    
    Returns
    -------
    tuple
        (lattice, positions, and numbers) used as inputs for spglib
        functions.
    """
    
    # Get lattice information
    lattice = system.box.vects
    
    # Get atomic positions (in reduced coordinates)
    positions = system.atoms_prop(key='pos', scale=True)
    
    # Get atom numbers (atypes)
    numbers = system.atoms.atype
    
    return (lattice, positions, numbers)