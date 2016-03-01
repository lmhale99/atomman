from DataModelDict import DataModelDict
import numpy as np
import atomman as am
import atomman.unitconvert as uc
import sys

def crystal(model):
    """Read in a data model containing a crystal-structure and return a System unit cell."""
    
    a_sys = DataModelDict(model).find('atom-system')

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
    sym_dict = {}
    scale = None
    
    prop = DataModelDict()
    
    if 'atom-sites' in a_sys:
        assert 'atom-properties' not in a_sys, 'atom-system model cannot have both atom-sites and atom-properties'
    
        for site in a_sys['atom-sites'].iterlist('site'):
            
            atype = site['type']

            if atype not in sym_dict:
                sym_dict[atype] = site.get('symbol', None)
            else:
                assert sym_dict[atype] == site.get('symbol', None), 'Inconsistent type-symbol pairings'
                
            #Add atype for site
            prop.append('atype', atype)
                        
            #read in pos for site and unit info
            pos = np.array(site['position']['value'], dtype='float64')
            unit = site['position'].get('unit', None)
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
                            
            #Add pos for site
            prop.append('pos', uc.set_in_units(pos, unit))
                
            #Add per-site properties       
            for property in site.iterlist('site-property'):
                prop.append(property['name'], uc.value_unit(property))
        
    elif 'atom-properties' in a_sys:
        a_props = a_sys['atom-properties']
        prop['atypes'] = uc.value_unit(a_props['type'])
        all_symbols = a_props['symbol'].list('value')
        for atype, sym in zip(prop['atypes'], all_symbols):
            if atype not in sym_dict:
                sym_dict[atype] = sym
            else:
                assert sym_dict[atype] == sym, 'Inconsistent type-symbol pairings'
        prop['pos'] = np.empty(len(prop['atypes']), 3)
        prop['pos'][:,0] = uc.value_unit(a_props['position-x'])
        prop['pos'][:,1] = uc.value_unit(a_props['position-y'])
        prop['pos'][:,2] = uc.value_unit(a_props['position-z'])
        
        unit = a_props['position-x'].get('unit', None)
        if unit == 'scaled':
            assert a_props['position-y'].get('unit', None) == 'scaled', 'Cannot mix box-scaled and absolute positions'
            assert a_props['position-z'].get('unit', None) == 'scaled', 'Cannot mix box-scaled and absolute positions'
            scale = True
        else:
            assert a_props['position-y'].get('unit', None) != 'scaled', 'Cannot mix box-scaled and absolute positions'
            assert a_props['position-z'].get('unit', None) != 'scaled', 'Cannot mix box-scaled and absolute positions'
            scale = False
        
        #Add additional properties       
        for property in site.iterlist('property'):
            prop[property['name']] = uc.value_unit(property)
        
    else:
        raise ValueError('atom-system model has neither atom-sites or atom-properties')
        
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