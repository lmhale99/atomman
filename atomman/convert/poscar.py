# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
import atomman.core.Atoms
import atomman.core.Box
import atomman.core.System
from ..compatibility import range
from ..tools import uber_open_rmode

def load(poscar, prop={}):
    """
    Reads a poscar-style coordination file for a system.
    
    Parameters
    ----------
    poscar : str or file-like object
        The POSCAR content to read
    prop : dict, optional
        Dictionary containing any extra per-atom properties to include.
    
    Returns
    -------
    system : atomman.System
        The system object associated with the poscar file.
    elements : list
        The list of elemental symbols corresponding to the system's atom types.
        Will be a list of None if symbol information is not in poscar.
    """
    
    # Read in all lines of the file
    with uber_open_rmode(poscar) as f:
        lines = f.readlines()
    
    # Interpret box information
    box_scale = float(lines[1])
    avect = np.array(lines[2].split(), dtype='float64') * box_scale
    bvect = np.array(lines[3].split(), dtype='float64') * box_scale
    cvect = np.array(lines[4].split(), dtype='float64') * box_scale
    box = atomman.core.Box(avect=avect, bvect=bvect, cvect=cvect)
    
    # Read in elements, number of types, and style info
    try:
        typenums = np.array(lines[5].split(), dtype='int64')
        elements = [None for n in range(len(typenums))]
        style = lines[6]
        start_i = 7
    except:
        elements = lines[5].split()
        typenums = np.array(lines[6].split(), dtype='int64')
        style = lines[7]
        start_i = 8
    
    for i in range(len(elements)):
        try:
            elements[i] = elements[i].decode('UTF-8')
        except:
            pass
    
    # Build atype list
    atype = np.array([], dtype='int64')
    for i in range(len(typenums)):
        atype = np.hstack((atype, np.full(typenums[i], i+1, dtype='int64')))
    
    # Check which coordinate style to use
    try:
        style = style.decode('UTF-8')
    except:
        pass
    if style[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    # Read in positions
    natoms = np.sum(typenums)
    pos = np.empty((natoms, 3), dtype='float64')
    count = 0
    for i in range(start_i, len(lines)):
        terms = lines[i].split()
        if len(terms) > 0:
            pos[count,:] = np.array(terms, dtype='float64')
        count += 1
    
    # Add atype, pos to prop dictionary
    prop['atype'] = atype
    prop['pos'] = pos
    
    atoms = atomman.core.Atoms(prop=prop)
    system = atomman.core.System(atoms=atoms, box=box, scale=scale)
    
    return system, elements

def dump(system, fname=None, header='', elements=None, style='direct',  box_scale=1.0, xf='%.13e'):
    """
    Generates a poscar-style coordination file for the system.
    
    Parameters
    ----------
    system : atomman.System
        The system whose coordinates you are saving
    fname : str, optional
        Path to file where POSCAR content is saved to.  If not given, a str of
        the content will be returned.
    header : str, optional
        The comment line to place at the top of the file. Default value is ''.
    elements : list of str, optional
        List of the element symbols that correspond to the atom types.  If not
        given, not element content will be included.
    style : str, optional 
        The poscar coordinate style.  Default value is 'direct'.
    box_scale : float, optional
        A universal scaling constant applied to the box vectors. Default value
        is 1.0.
     xf : str, optional
        c-style format for printing the floating point numbers. Default value
        is '%.13e'.
        
    Returns
    -------
    poscar_str : str
        String of the poscar object (only returned if fname is not given).
    """
    assert '\n' not in header, 'header can only be one line'
    assert '\n' not in style, 'style can only be one line'
    
    threexf = xf + ' ' + xf + ' ' + xf
    
    # Scale box vectors and write out the values
    vects = system.box.vects / box_scale
    poscar_string = '\n'.join([header,
                               xf % box_scale,
                               threexf % tuple(vects[0]),
                               threexf % tuple(vects[1]),
                               threexf % tuple(vects[2])])
    
    # Write element tags if they are given
    if elements is not None:
        if not isinstance(elements, (list, tuple)):
            elements = [elements]
        assert len(elements) == system.natypes, 'length of elements differs from number of atom types'
        poscar_string += '\n' + ' '.join(elements)
    
    # Count how many atoms of each type
    atype = system.atoms.atype
    poscar_string += '\n' 
    uatype, counts = np.unique(atype, return_counts=True)

    for i in range(1, int(uatype.max()+1)):
        count = counts[uatype==i]
        if count == []:
            count = 0
        else:
            count = count[0]
        poscar_string += '%i ' % count
        
    # Check which coordinate style to use
    poscar_string += '\n' + style
    if style[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    # Write out positions
    pos = system.atoms_prop(key='pos', scale=scale)
    for a in range(1, system.natypes+1):
        for p in pos[atype==a]:
            poscar_string += '\n'+ threexf % tuple(p)
    
    # Save to the file name
    if fname is not None:
        with open(fname, 'w') as f:
            f.write(poscar_string)
    
    # Return as a string
    else:
        return poscar_string
   