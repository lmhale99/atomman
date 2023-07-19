# coding: utf-8
# Standard Python libraries
import io
from math import ceil
from typing import Union

def dump(system,
         f: Union[str, io.IOBase, None] = None,
         ignoresymbols: bool = False):
    """
    Dumps the atomic configuration to the pdb format.  Since atomman only tracks
    atomic positions, this simply sets the atomic configurations and cell
    dimensions. Also, note that coordinates will be shifted such that the cell
    origin is at (0,0,0).

    Parameters
    ----------
    system : atomman.System
        The system whose coordinates you are saving.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    ignoresymbols : bool, optional
        Setting this to True will save all atom types as unknowns rather than
        as their element symbols.  This is useful in tricking some pdb-based
        tools into recognizing the atoms if the model symbols differ from
        elemental symbols or the element is H. 
    """
    # Build CRYST1 from box information
    a = system.box.a
    b = system.box.b
    c = system.box.c
    alpha = system.box.alpha
    beta = system.box.beta
    gamma = system.box.gamma
    lines = [f'CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1']

    # Extract atomic information
    atype = system.atoms.atype - 1
    pos = system.atoms.pos - system.box.origin
    natoms = system.natoms
    symbols = system.symbols

    nmodels = ceil(system.natoms / 99999)
    for model in range(nmodels):

        lines.append(f'MODEL     {model+1}')

        start = model * 99999
        end = (model + 1) * 99999
        if end > natoms:
            end = natoms
        nitems = end - start

        for atom in range(nitems):
            id = atom + 1
            aid = start + atom
            name = symbols[atype[aid]]
            if ignoresymbols or name is None:
                name = 'UNX'
                symbol = f'X{atype[aid]+1}'
            else:
                symbol = name.upper()
            lines.append(f'HETATM{id:5} {name:4} MOL     1    {pos[aid,0]:8.3f}{pos[aid,1]:8.3f}{pos[aid,2]:8.3f}  1.00  0.00          {symbol:>2}')
        lines.append('ENDMDL')
    lines.append('')

    pdb_string = '\n'.join(lines)

    # Save to the file-like object
    if hasattr(f, 'write'):
        f.write(pdb_string)

    # Save to the file name
    elif f is not None:
        with open(f, 'w', encoding='UTF-8') as fp:
            fp.write(pdb_string)

    # Return as a string
    else:
        return pdb_string
