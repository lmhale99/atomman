from DataModelDict import DataModelDict
import numpy as np
import atomman as am
import atomman.unitconvert as uc
import sys

def crystal(model):
    """Read in a data model containing a crystal-structure and return a System unit cell."""
    
    cell = DataModelDict(model).find('cell')
    
    #identify the crystal system
    c_system = cell.keys()[0]
    cell = cell[c_system]
    
    if c_system == 'cubic': 
        a = b = c = uc.value_unit(cell['a'])
        alpha = beta = gamma = 90.0
        
    elif c_system == 'hexagonal':
        a = b = uc.value_unit(cell['a'])
        c =     uc.value_unit(cell['c'])
        alpha = beta = 90.0
        gamma =       120.0
        
    elif c_system == 'tetragonal':
        a = b = uc.value_unit(cell['a'])
        c =     uc.value_unit(cell['c'])
        alpha = beta = gamma = 90.0
        
    elif c_system == 'trigonal' or c_system == 'rhombohedral':
        a = b = c = uc.value_unit(cell['a'])
        alpha = beta = gamma = cell['alpha']
        
    elif c_system == 'orthorhombic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = beta = gamma = 90.0

    elif c_system == 'monoclinic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = gamma = 90.0
        beta = cell['beta']
        
    elif c_system == 'triclinic':
        a = uc.value_unit(cell['a'])
        b = uc.value_unit(cell['b'])
        c = uc.value_unit(cell['c'])
        alpha = cell['alpha']
        beta =  cell['beta']
        gamma = cell['gamma']
        
    box = am.Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    
    #create list of atoms and list of elements
    atoms = []
    sym_dict = {}
    scale = None
    
    prop = DataModelDict()
    for atom in cell.iterlist('atom'):
        
        atype = atom['constituent']['component']

        if atype not in sym_dict:
            sym_dict[atype] = atom['constituent'].get('symbol', None)
        else:
            assert sym_dict[atype] == atom['constituent'].get('symbol', None), 'Inconsistent constituent component-symbols'
            
        #Add atype for atom
        prop.append('atype', atype)
                    
        #read in pos for atom and unit info
        pos = np.array(atom['position']['value'], dtype='float64')
        unit = atom['position'].get('unit', None)
        if unit == 'scaled':
            unit = None
            if scale:
                pass
            elif scale is None:
                scale = True
            else:
                raise ValueError('Cannot mix box-scaled and absolute positions')
        else:
            scale = False
                        
        #Add pos for atom
        prop.append('pos', uc.set_in_units(pos, unit))
            
        #Add per-atom properties       
        for property in atom.iterlist('atom-property'):
            prop.append(property['name'], uc.value_unit(property))
    
        atype += 1
    atoms = am.Atoms(natoms=len(prop['atype']), prop=prop)
    system = am.System(box=box, atoms=atoms, scale=scale)
    
    symbols = [None for i in xrange(max(sym_dict.keys()))]
    for i in sym_dict.keys():
        symbols[i-1] = sym_dict[i]
    
    return system, symbols        
    
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            ucell, symbols = crystal(f)
        print 'symbols =', symbols
        print ucell