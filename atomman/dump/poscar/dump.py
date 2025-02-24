# coding: utf-8

# http://www.numpy.org/
import io
from typing import Optional, Union
import numpy as np

def dump(system,
         f: Union[str, io.IOBase, None] = None,
         header: str = '',
         symbols: Optional[tuple] = None,
         coordstyle: str = 'direct',
         box_scale: float = 1.0,
         float_format: str = '%.13e') -> Optional[str]:
    """
    Generates a poscar-style coordination file for the system.
    
    Parameters
    ----------
    system : atomman.System
        The system whose coordinates you are saving
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    header : str, optional
        The comment line to place at the top of the file. Default value is ''.
    symbols : tuple, optional
        List of the element symbols that correspond to the atom types.  If not
        given, will use system.symbols if set, otherwise no element content
        will be included.
    coordstyle : str, optional
        The poscar coordinate style to use: 'cartesian' or 'direct' (i.e.
        box relative).  Default value is 'direct'.
    box_scale : float, optional
        A universal scaling constant applied to the box vectors. Default value
        is 1.0.
     float_format : str, optional
        c-style format for printing the floating point numbers. Default value
        is '%.13e'.
        
    Returns
    -------
    poscar_str : str
        String of the poscar object (only returned if fname is not given).
    """
    assert '\n' not in header, 'header can only be one line'
    assert '\n' not in coordstyle, 'coordstyle can only be one line'
    
    threexf = float_format + ' ' + float_format + ' ' + float_format
    
    # Scale box vectors and write out the values
    vects = system.box.vects / box_scale
    poscar_string = '\n'.join([header,
                               float_format % box_scale,
                               threexf % tuple(vects[0]),
                               threexf % tuple(vects[1]),
                               threexf % tuple(vects[2])])
    
    # Use system.symbols if set
    if symbols is None:
        if None not in system.symbols:
            symbols = system.symbols
    
    # Write symbols tags if they are given
    if symbols is not None:
        if not isinstance(symbols, (list, tuple)):
            symbols = [symbols]
        if len(symbols) != system.natypes:
            raise ValueError('length of symbols differs from number of atom types')
        poscar_string += '\n' + ' '.join(symbols)
    
    # Count how many atoms of each type
    atype = system.atoms.atype
    poscar_string += '\n' 
    uatype, counts = np.unique(atype, return_counts=True)
    
    for i in range(1, int(uatype.max()+1)):
        count = counts[uatype==i]
        if count.size == 0:
            count = 0
        else:
            count = count[0]
        poscar_string += '%i ' % count
        
    # Check which coordinate style to use
    poscar_string += '\n' + coordstyle
    if coordstyle[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    # Write out positions
    pos = system.atoms_prop(key='pos', scale=scale)
    for a in range(1, system.natypes+1):
        for p in pos[atype==a]:
            poscar_string += '\n'+ threexf % tuple(p)
    
    # Save to the file-like object
    if hasattr(f, 'write'):
        f.write(poscar_string)
    
    # Save to the file name
    elif f is not None:
        with open(f, 'w', encoding='UTF-8') as fp:
            fp.write(poscar_string)
    
    # Return as a string
    else:
        return poscar_string