from DataModelDict import DataModelDict
import numpy as np
import atomman as am
import atomman.unitconvert as uc
import sys

def crystal(model):
    """Read in a data model containing a crystal-structure and return a System unit cell."""
    
    struct = DataModelDict(model).find('crystal-structure')
    
    #identify the crystal system
    c_system = struct['cell'].keys()[0]
    cell = struct.find(c_system)
    
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
    symbols = []
    atype = 1
    scale = None
    
    prop = DataModelDict()
    for site in struct.iterlist('site'):
        atype_props = {}
        
        #set atom type values
        assert atype == site['constituent']['component'], 'invalid site component'
        try:
            symbol = site['constituent']['symbol']
        except:
            try:
                symbol = site['constituent']['element']
            except:
                symbol = None
        symbols.append(symbol)
        
        for property in site['constituent'].iterlist('atom-property'):
            atype_props[property['name']] = uc.value_unit(property)

        #set per-atom values
        for atom in site.iterlist('atom'):
            
            #Add atype for atom
            prop.append('atype', atype)
                    
            #read in pos for atom and unit info
            pos = np.array(atom['equivalent-position']['value'], dtype='float64')
            unit = atom['equivalent-position'].get('unit', None)
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
            
            #Add atom-type properties
            for k, v in atype_props.iteritems():
                prop.append(k, v)
                
            #Add per-atom properties       
            for property in site['constituent'].iterlist('atom-property'):
                prop.append(property['name'], uc.value_unit(property))
        
        atype += 1
    atoms = am.Atoms(natoms=len(prop['atype']), prop=prop)
    system = am.System(box=box, atoms=atoms, scale=scale)
    
    return system, symbols        
    
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            ucell = crystal(f)
        print 'symbols =', ucell[1]
        print ucell[0]