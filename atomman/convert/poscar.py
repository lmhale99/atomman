import atomman as am
import numpy as np

from atomman.tools import uber_open_rmode

def load(poscar):
    """
    Reads a poscar-style coordination file for a system.
    Returns an atomman.System, and a list of elements if the file gives them. 
    """
    #Read in all lines of the file
    with uber_open_rmode(poscar) as f:
        lines = f.read().split('\n')
    
    #Interpret box information
    box_scale = float(lines[1])
    avect = np.array(lines[2].split(), dtype='float64') * box_scale
    bvect = np.array(lines[3].split(), dtype='float64') * box_scale
    cvect = np.array(lines[4].split(), dtype='float64') * box_scale
    box = am.Box(avect=avect, bvect=bvect, cvect=cvect)
    
    #Read in elements, number of types, and style info
    try:
        typenums = np.array(lines[5].split(), dtype='int32')
        elements = [None for n in xrange(len(typenums))]
        style = lines[6]
        start_i = 7
    except:
        elements = lines[5].split()
        typenums = np.array(lines[6].split(), dtype='int32')
        style = lines[7]
        start_i = 8

    #Build atype list
    atype = np.array([], dtype='int32')
    for i in xrange(len(typenums)):
        atype = np.hstack((atype, np.full(typenums[i], i+1, dtype='int32')))        
    
    #Check which coordinate style to use
    if style[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    #Read in positions
    natoms = np.sum(typenums)
    pos = np.empty((natoms, 3), dtype='float64')
    count = 0
    for i in xrange(start_i, len(lines)):
        terms = lines[i].split()
        if len(terms) > 0:
            pos[count,:] = np.array(terms, dtype='float64')
        count += 1 
    
    atoms = am.Atoms(natoms=natoms, prop={'atype':atype, 'pos':pos})
    system = am.System(atoms=atoms, box=box, scale=scale)
    
    return system, elements
    
        

def dump(system, fp=None, fname=None, header='', elements=None, style='direct',  box_scale=1.0):
    """
    Generates a poscar-style coordination file for the system.
    
    Arguments:
    fname -- name of the file to save to.
    system -- the atomman.System whose coordinates you are saving.
    
    Optional Keyword Arguments:
    header -- string of the comment line to place at the top of the file. Default is ''.
    elements -- list of the elements that correspond to the atom types.
    style -- the poscar coordinate style: Direct or Cartesian.
    box_scale -- a universal scaling constant applied to the box vectors. Default is 1.0.
    """
    assert '\n' not in header, 'header can only be one line'
    assert '\n' not in style, 'style can only be one line'
    assert fname is None or fp is None, 'fp and fname cannot both be given'
    
    box_scale = float(box_scale)
    
    #scale box vectors and write out the values
    vects = system.box.vects / box_scale
    poscar_string = '\n'.join([header,
                               repr(box_scale),
                               '%s %s %s' % (repr(vects[0,0]), repr(vects[0,1]), repr(vects[0,2])),
                               '%s %s %s' % (repr(vects[1,0]), repr(vects[1,1]), repr(vects[1,2])),
                               '%s %s %s' % (repr(vects[2,0]), repr(vects[2,1]), repr(vects[2,2]))])
    
    #Write element tags if they are given
    if elements is not None:
        if not isinstance(elements, (list, tuple)):  elements = [elements]
        assert len(elements) == system.natypes, 'length of elements differs from number of atom types'
        poscar_string += '\n' + ' '.join(elements)
    
    #Count how many atoms of each type
    atype = system.atoms_prop(key = 'atype')
    poscar_string += '\n' 
    for bc in np.bincount(atype)[1:]:
        poscar_string += '%i ' % bc
        
    #Check which coordinate style to use
    poscar_string += '\n' + style
    if style[0] in 'cCkK':
        scale = False
    else:
        scale = True
    
    #Write out positions
    pos = system.atoms_prop(key='pos', scale=scale)
    for a in xrange(1, system.natypes+1):
        for p in pos[atype==a]:
            poscar_string += '\n%s %s %s' % (repr(p[0]), repr(p[1]), repr(p[2]))
    
    #Save to the file pointer
    if fp is not None:
        fp.write(poscar_string) 
    
    #Save to the file name
    elif fname is not None:
        with open(fname, 'w') as f:
            f.write(poscar_string)
    
    #Return as a string
    else:
        return poscar_string
   