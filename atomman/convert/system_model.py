from DataModelDict import DataModelDict
import numpy as np
import atomman as am
import atomman.unitconvert as uc
import sys
import os

def load(model, key='atomic-system', index=0):
    """Read in a data model containing a crystal-structure and return a System unit cell."""
    
    if isinstance(model, (str, unicode)) and os.path.isfile(model):
        with open(model) as f:
            model = f.read()
    
    #Pull system model out of data model using key and index
    a_sys = DataModelDict(model).finds(key)
    if len(a_sys) == 0:
        raise KeyError(key + ' not found in model')
    try:
        a_sys = a_sys[index]
    except:
        raise IndexError('Invalid index ' + str(index) + ' for model key ' + key)

    #identify the crystal system
    c_system = a_sys['cell'].keys()[0]
    cell = a_sys['cell'][c_system]
    
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
    scale = None
    
    prop = DataModelDict()
    
    all_atypes = np.array(a_sys.finds('component'))
    all_symbols = np.array(a_sys.finds('symbol'))
    all_elements = np.array(a_sys.finds('element'))
    
    if len(all_atypes) == 0:
        if len(all_symbols) != 0:
            symbols, atypes = np.unique(all_symbols, return_inverse)
        elif len(all_elements) != 0:
            symbols, atypes = np.unique(all_elements, return_inverse)
        else:
            raise ValueError('No atom components, symbols or elements listed')
            
    else:
        atypes = all_atypes
        symbols = [None for i in xrange(max(all_atypes))]
        
        if len(all_elements) != 0 and len(all_symbols) == 0:
            all_symbols = all_elements
        
        if len(all_symbols) != 0:
            assert len(all_symbols) == len(atypes)
            sym_dict = {}
            for atype, symbol in zip(atypes, all_symbols):
                if atype not in sym_dict:
                    sym_dict[atype] = symbol
                else:
                    assert sym_dict[atype] == symbol
            
            for atype, symbol in sym_dict.iteritems():
                symbols[atype-1] = symbol
    
    prop['atype'] = atypes
    prop['pos'] = np.zeros((len(prop['atype']), 3), dtype='float64')
    count = 0
    
    pos_units = []
    for atom in a_sys.iteraslist('atom'):
                      
        #read in pos for atom and unit info
        prop['pos'][count] = uc.value_unit(atom['position'])
        pos_units.append(atom['position'].get('unit', None))
            
        #Add per-atom properties       
        for property in atom.iteraslist('property'):
            if property['name'] not in prop:
                value = uc.value_unit(property)
                prop[property['name']] = np.zeros((len(prop['atype']), len(value)), dtype=value.dtype)
                
            prop[property['name']][count] = uc.value_unit(property)
        count += 1
    
    pos_unit = np.unique(pos_units)
    assert len(pos_unit) == 1, 'Mixed units for positions'
    if pos_unit[0] == 'scaled': scale=True
    else:                       scale=False
    
    atoms = am.Atoms(natoms=len(prop['atype']), prop=prop)
    system = am.System(box=box, atoms=atoms, scale=scale)
    
    return system, symbols        
    
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            ucell, symbols = crystal(f)
        print 'symbols =', symbols
        print ucell