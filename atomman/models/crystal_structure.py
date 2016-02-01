from DataModelDict import DataModelDict 
import numpy as np
import atomman

def crystal_structure(model, **kwargs):
    """Build System based on data model containing crystal-structure info.
    
    Argument:
    model -- string or file-like object json/xml containing a single crystal-structure branch 
    
    Keyword Arguments:
    a, b, c, alpha, beta, gamma -- optional lattice parameters to overwrite what is in the datamodel.
    
    Note that the structure's crystal system influences which lattice parameter arguments can be specified, and whether other lattice parameters are scaled."""
    
    #extract crystal-structure information from the data model
    crystal = DataModelDict(model).find('crystal-structure')
    assert len(crystal) == 1, 'Exactly one crystal structure must be listed in data model'
    crystal = crystal[0]
    
    #identify the crystal_system
    crystal_system = None
    for sys_type in ('cubic', 'tetragonal', 'hexagonal', 'orthorhombic', 'rhombohedral', 'monoclinic', 'trigonal', 'triclinic'):
        try:
            cell = crystal['cell'][sys_type]
            crystal_system = sys_type
            break
        except:
            pass
    assert crystal_system is not None, 'No crystal system found!'
    
    #extract default cell lengths
    default = {}
    default['a'] = atomman.value_unit(cell['a']) 
    if crystal_system in ('orthorhombic', 'monoclinic', 'triclinic'):
        default['b'] = atomman.value_unit(cell['b']) 
        default['c'] = atomman.value_unit(cell['c']) 
    elif crystal_system in ('tetragonal', 'hexagonal'):
        default['c'] = atomman.value_unit(cell['c']) 
    
    #extract default cell angles
    if crystal_system in ('trigonal', 'rhomohedral'):
        default['alpha'] = cell['alpha']['value']
    elif crystal_system == 'monoclinic':
        default['beta'] =  cell['beta']['value']
    elif crystal_system == 'triclinic':
        default['alpha'] = cell['alpha']['value']
        default['beta'] = cell['beta']['value']
        default['gamma'] = cell['gamma']['value']
    
    #check cell angles units
    try:    a_unit = cell['alpha']['unit']
    except: a_unit = None
    try:    b_unit = cell['beta']['unit']
    except: b_unit = None
    try:    c_unit = cell['gamma']['unit']
    except: c_unit = None
    assert a_unit == 'degree' or a_unit is None, 'cell angles must be given in degrees'
    assert b_unit == 'degree' or b_unit is None, 'cell angles must be given in degrees'
    assert c_unit == 'degree' or c_unit is None, 'cell angles must be given in degrees'
    
    #make a box based on default and argument cell parameters
    box = make_box(crystal_system, default, kwargs) 
    
    #create list of atoms and list of elements
    atoms = []
    el_list = []
    atype = 1
    scale = None
    
    sites = crystal['site']
    if not isinstance(sites, (list)):
        sites = [sites]
    
    prop = DataModelDict()
    for site in sites:
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
        try:
            properties = site['constituent']['atom-property']
            if not isinstance(properties, list):
                properties = [properties]
        except:
            properties = []       
        for property in properties:
            atype_props[property['name']] = atomman.value_unit(property)
        el_list.append(symbol)
        
        #set per-atom values
        if not isinstance(site['atom'], list):
            site['atom'] = [site['atom']]
        for atom in site['atom']:
            
            #Add atype for atom
            try:    prop['atype'].append(atype)
            except: prop['atype'] = [atype]
        
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
            pos = atomman.set_in_units(pos, unit)
            
            #Add pos for atom
            try:    prop['pos'].append(pos)
            except: prop['pos'] = [pos]
            
            #Add atom-type properties
            for k, v in atype_props.iteritems():
                try:    prop[k].append(v)
                except: prop[k] = [v]
                
            #Add per-atom properties
            try:
                properties = site['constituent']['atom-property']
                if not isinstance(properties, list):
                    properties = [properties]
            except:
                properties = []       
            for property in properties:
                try:    prop[property['name']].append(atomman.value_unit(property))
                except: prop[property['name']] = [atomman.value_unit(property)]
        
        atype += 1
    atoms = atomman.Atoms(natoms=len(prop['atype']), prop=prop)
    system = atomman.System(box=box, atoms=atoms, scale=scale)
    
    return system, el_list        
    
def make_box(crystal_system, default, arg):
    
    if crystal_system == 'cubic':
        if arg['a'] is None:
            assert arg['b'] is None and arg['c'] is None, 'define lattice parameters in order'
            a = b = c = default['a']
        else:
            assert (arg['a'] == arg['b'] or arg['b'] is None) and (arg['a'] == arg['c'] or arg['c'] is None), 'cubic lengths must match'
            a = b = c = arg['a']
            
        assert arg['alpha'] is None or arg['alpha'] == 90.0, 'cubic angles must be 90'
        assert arg['beta']  is None or arg['beta']  == 90.0, 'cubic angles must be 90'
        assert arg['gamma'] is None or arg['gamma'] == 90.0, 'cubic angles must be 90'
        alpha = beta = gamma = 90.0
            
    elif crystal_system == 'tetragonal':
        if arg['a'] is None:
            assert arg['b'] is None, 'define lattice parameters in order'
            a = b = default['a']
        else:
            assert arg['b'] is None or arg['a'] == arg['b'], 'tetragonal lengths a and b must match'
            a = b = arg['a']
        if arg['c'] is None:
            c = default['c'] * a / default['a']
        else:
            c = arg['c']    
            
        assert arg['alpha'] is None or arg['alpha'] == 90.0, 'tetragonal angles must be 90'
        assert arg['beta']  is None or arg['beta']  == 90.0, 'tetragonal angles must be 90'
        assert arg['gamma'] is None or arg['gamma'] == 90.0, 'tetragonal angles must be 90'
        alpha = beta = gamma = 90.0
        
    elif crystal_system == 'hexagonal':
        if arg['a'] is None:
            assert arg['b'] is None, 'define lattice parameters in order'
            a = b = default['a']
        else:
            assert arg['b'] is None or arg['a'] == arg['b'], 'hexagonal lengths a and b must match'
            a = b = arg['a']
        if arg['c'] is None:
            c = default['c'] * a / default['a']
        else:
            c = arg['c']  
            
        assert arg['alpha'] is None or arg['alpha'] == 90.0, 'hexagonal angles must be 90'
        assert arg['beta']  is None or arg['beta']  == 90.0, 'hexagonal angles must be 90'
        assert arg['gamma'] is None or arg['gamma'] == 120.0,'hexagonal angles must be 90'
        alpha = beta = 90.0
        gamma = 120.0
    
    elif crystal_system == 'orthorhombic':
        if arg['a'] is None:
            a = default['a']
        else:
            a = arg['a']
        if arg['b'] is None:
            b = default['b'] * a / default['a']
        else:
            b = arg['b']     
        if arg['c'] is None:
            c = default['c'] * a / default['a']
        else:
            c = arg['c'] 
        
        assert arg['alpha'] is None or arg['alpha'] == 90.0, 'orthorhombic angles must be 90'
        assert arg['beta']  is None or arg['beta']  == 90.0, 'orthorhombic angles must be 90'
        assert arg['gamma'] is None or arg['gamma'] == 90.0, 'orthorhombic angles must be 90'
        alpha = beta = gamma = 90.0
    
    elif crystal_system == 'rhombohedral' or crystal_system == 'trigonal':
        if arg['a'] is None:
            assert arg['b'] is None and arg['c'] is None, 'define lattice parameters in order'                    
            a = b = c = default['a']
        else:
            assert (arg['a'] == arg['b'] or arg['b'] is None) and (arg['a'] == arg['c'] or arg['c'] is None), 'rhombohedral lengths must match'
            a = b = c = arg['a']
        
        if arg['alpha'] is None:
            assert arg['beta'] is None and arg['gamma'] is None, 'define lattice angles in order'                    
            alpha = beta = gamma = default['alpha']
        else:
            assert (arg['alpha'] == arg['beta'] or arg['beta'] is None) and (arg['alpha'] == arg['gamma'] or arg['gamma'] is None), 'rhombohedral angles must match'
            alpha = beta = gamma = arg['alpha']
        assert alpha < 120, 'Non-standard rhombohedral angle'
    
    elif crystal_system == 'monoclinic':
        if arg['a'] is None:
            a = default['a']
        else:
            a = arg['a']
        if arg['b'] is None:
            b = default['b'] * a / default['a']
        else:
            b = arg['b']     
        if arg['c'] is None:
            c = default['c'] * a / default['a']
        else:
            c = arg['c'] 
        
        assert arg['alpha'] is None or arg['alpha'] == 90.0, 'orthorhombic angles must be 90'
        assert arg['gamma'] is None or arg['gamma'] == 90.0, 'orthorhombic angles must be 90'
        alpha = gamma = 90.0
        if arg['beta'] is None:
            beta = default['beta']
        else:
            beta = arg['beta']
        assert beta > 90, 'Non-standard monoclinic angle'

    elif crystal_system == 'triclinic':
        if arg['a'] is None:
            a = default['a']
        else:
            a = arg['a']
        if arg['b'] is None:
            b = default['b'] * a / default['a']
        else:
            b = arg['b']     
        if arg['c'] is None:
            c = default['c'] * a / default['a']
        else:
            c = arg['c'] 

        if arg['alpha'] is None:
            alpha = default['alpha']
        else:
            alpha = arg['alpha']
        if arg['beta'] is None:
            beta = default['beta']
        else:
            beta = arg['beta']     
        if arg['gamma'] is None:
            gamma = default['gamma']
        else:
            gamma = arg['gamma'] 

    return atomman.Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)        
        
        
        
        
        
        