import atomman as am
import atomman.unitconvert as uc
from DataModelDict import DataModelDict
import numpy as np
from copy import deepcopy

def load(fname, prop_info=None):
    """
    Read a LAMMPS-style dump file and return a System.
    
    Argument:
    fname -- name (and location) of file to read data from.
    
    Keyword Argument:
    prop_info -- DataModelDict for relating the per-atom properties to/from the dump file and the System. Will create a default json instance <fname>.json if prop_info is not given and <fname>.json doesn't already exist.
    """
    
    #read in prop_info if supplied
    if prop_info is not None:
        prop_info = DataModelDict(prop_info)
    
    #check for default prop_info file
    else:
        try:
            with open(fname+'.json') as fj:
                prop_info = DataModelDict(fj)
        except:
            prop_info = None
            box_unit = None
   
    #read box_unit if specified in prop_info
    if prop_info is not None:
        prop_info = prop_info.find('LAMMPS-dump-atoms_prop-relate')
        box_unit = prop_info['box_prop'].get('unit', None)

    with open(fname, 'r') as f:
        pbc = None
        box = None
        natoms = None
        
        readnatoms = False
        readatoms = False
        readtimestep = False
        acount = 0
        bcount = 3
        
        #loop over all lines in file
        for line in f:
            terms = line.split()
            if len(terms) > 0:
    
                #read atomic values if time to do so
                if readatoms:
                    #sort values by a_id and save to prop_vals
                    a_id = long(terms[id_index]) - 1
                    prop_vals[a_id] = terms
                    acount += 1
                    
                    #save values to sys once all atoms read in
                    if acount == natoms:
                        readatoms = False
                        
                        #cycle over the defined atoms_prop in prop_info
                        for prop, p_keys in prop_info['atoms_prop'].iteritems():
                            #set default keys
                            dtype = p_keys.get('dtype', None)
                            shape = p_keys.get('shape', None)
                            shape = (natoms,) + np.empty(shape).shape
                            
                            value = np.empty(shape)
                            
                            #cycle over the defined LAMMPS-attributes in prop_info
                            for attr, a_keys in prop_info['LAMMPS-attribute'].iteritems():
                                
                                #cycle over list of relations for each LAMMPS-attribute
                                for relation in a_keys.iterlist('relation'):
                                    
                                    #if atoms prop and relation prop match
                                    if relation['prop'] == prop:
                                        #get unit and scale info
                                        unit = relation.get('unit', None)

                                        if unit == 'scaled':
                                            unit = None
                                            scale = True
                                        else:
                                            scale = False
                                        
                                        #find index of attribute in name_list
                                        a_index = name_list.index(attr) 
                                        #check if relation has index listed
                                        try:
                                            index = relation['index']
                                            if isinstance(index, list):
                                                index = (Ellipsis,) + tuple(index)
                                            else:
                                                index = (Ellipsis,) + (index,)
                                                
                                            value[index] = prop_vals[:, a_index]
                                        #scalar if no index
                                        except:
                                            value[:] = prop_vals[:, a_index]
                            #test if values are ints if dtype not specified
                            if dtype is None and np.allclose(np.asarray(value, dtype=int), value):
                                value = np.asarray(value, dtype=int)
                            else:
                                value = np.asarray(value, dtype=dtype)
                            
                            #save prop values to system
                            system.atoms_prop(key=prop, value=uc.set_in_units(value, unit), scale=scale)
                
                #read number of atoms if time to do so
                elif readnatoms:                
                    natoms = int(terms[0])
                    readnatoms = False
                
                elif readtimestep:
                    timestep = int(terms[0])
                    readtimestep = False
                
                #read x boundary condition values if time to do so
                elif bcount == 0:
                    xlo = uc.set_in_units(float(terms[0]), box_unit)
                    xhi = uc.set_in_units(float(terms[1]), box_unit)
                    if len(terms) == 3:
                        xy = uc.set_in_units(float(terms[2]), box_unit)
                    bcount += 1
                    
                #read y boundary condition values if time to do so
                elif bcount == 1:
                    ylo = uc.set_in_units(float(terms[0]), box_unit)
                    yhi = uc.set_in_units(float(terms[1]), box_unit)
                    if len(terms) == 3:
                        xz = uc.set_in_units(float(terms[2]), box_unit)
                    bcount += 1
                    
                #read z boundary condition values if time to do so
                elif bcount == 2:
                    zlo = uc.set_in_units(float(terms[0]), box_unit)
                    zhi = uc.set_in_units(float(terms[1]), box_unit)
                    if len(terms) == 3:
                        yz = uc.set_in_units(float(terms[2]), box_unit)
                        xlo = xlo - min((0.0, xy, xz, xy + xz))
                        xhi = xhi - max((0.0, xy, xz, xy + xz))
                        ylo = ylo - min((0.0, yz))
                        yhi = yhi - max((0.0, yz))
                        box = am.Box(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi, xy=xy, xz=xz, yz=yz)
                    else:
                        box = am.Box(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi) 
                    bcount += 1
                
                #if not time to read value, check the ITEM: header information
                else:
                    
                    #only consider ITEM: lines
                    if terms[0] == 'ITEM:':
                    
                        #ITEM: TIMESTEP indicates it is time to read the timestep
                        if terms[1] == 'TIMESTEP':
                            readtimestep = True                        
                        
                        #ITEM: NUMBER indicates it is time to read natoms
                        elif terms[1] == 'NUMBER':
                            readnatoms = True
                            
                        #ITEM: BOX gives pbc and indicates it is time to read box parameters
                        elif terms[1] == 'BOX':
                            pbc = [True, True, True]
                            for i in xrange(3):
                                if terms[i + len(terms) - 3] != 'pp':
                                    pbc[i] = False
                            bcount = 0
                            
                        #ITEM: ATOMS gives list of per-Atom property names and indicates it is time to read atomic values      
                        elif terms[1] == 'ATOMS':
                            assert box is not None,  'Box information not found'
                            assert natoms is not None,  'Number of atoms not found'
                            
                            #read list of property names
                            name_list = terms[2:]
                            id_index = name_list.index('id')
                            
                            #create empty array for reading property values
                            prop_vals = np.empty((natoms, len(name_list)))
                            
                            #create and save default prop_info Data Model if needed
                            if prop_info is None:                                                       
                                prop_info = __prop_info_default_load(name_list)
                                with open(fname+'.json', 'w') as fj:
                                    prop_info.json(fp=fj, indent=4)
                                prop_info = prop_info.find('LAMMPS-dump-atoms_prop-relate')
                            
                            #create system and flag that it is time to read data
                            system = am.System(atoms=am.Atoms(natoms=natoms), box=box, pbc=pbc)
                            system.prop['timestep'] = timestep
                            readatoms = True 
    return system      

def dump(fname, system, prop_info=None, xf='%.13e'):
    """
    Write a LAMMPS-style dump file from a System.
    
    Arguments:
    fname -- name (and location) of file to save data to.
    system -- System to write to the dump file.
    
    Keyword Arguments:
    prop_info -- DataModelDict for relating the per-atom properties to/from the dump file and the System. Will create a default json instance <fname>.json if prop_info is not given and <fname>.json doesn't already exist.
    xf -- c-style format for printing the floating point numbers. Default is '%.13e'.
    """  

    #create or read prop_info Data Model
    if prop_info is None:
        try:
            with open(fname+'.json') as fj:
                prop_info = DataModelDict(fj)
        except:                            
            prop_info = __prop_info_default_dump(system)
            with open(fname+'.json', 'w') as fj:
                prop_info.json(fp=fj, indent=4)
    else:
        prop_info = DataModelDict(prop_info)
    
    #read box_unit if specified in prop_info
    prop_info = prop_info.find('LAMMPS-dump-atoms_prop-relate')
    box_unit = prop_info['box_prop'].get('unit', None)
        
    #open fname
    with open(fname, 'w') as f:

        #write timestep info
        f.write('ITEM: TIMESTEP\n')
        try:
            f.write('%i\n'%system.prop['timestep'])
        except:
            f.write('0\n')
        
        #write number of atoms
        f.write('ITEM: NUMBER OF ATOMS\n')
        f.write('%i\n' % ( system.natoms ))
        
        #write system boundary info for an orthogonal box
        if system.box.xy == 0.0 and system.box.xz == 0.0 and system.box.yz == 0.0:
            f.write('ITEM: BOX BOUNDS')
            for i in xrange(3):
                if system.pbc[i]:
                    f.write(' pp')
                else:
                    f.write(' fm')
            f.write('\n')
            
            f.write('%f %f\n' % ( uc.get_in_units(system.box.xlo, box_unit), uc.get_in_units(system.box.xhi, box_unit) )) 
            f.write('%f %f\n' % ( uc.get_in_units(system.box.ylo, box_unit), uc.get_in_units(system.box.yhi, box_unit) ))
            f.write('%f %f\n' % ( uc.get_in_units(system.box.zlo, box_unit), uc.get_in_units(system.box.zhi, box_unit) ))
        
        #write system boundary info for a triclinic box
        else:
            f.write('ITEM: BOX BOUNDS xy xz yz')
            for i in xrange(3):
                if system.pbc(i):
                    f.write(' pp')
                else:
                    f.write(' fm')
            f.write('\n')

            xlo_bound = uc.get_in_units(system.box.xlo, box_unit) + uc.get_in_units(min(( 0.0, system.box.xy, system.box.xz, system.box.xy + system.box.xz)), box_unit)
            xhi_bound = uc.get_in_units(system.box.xhi, box_unit) + uc.get_in_units(max(( 0.0, system.box.xy, system.box.xz, system.box.xy + system.box.xz)), box_unit)
            ylo_bound = uc.get_in_units(system.box.ylo, box_unit) + uc.get_in_units(min(( 0.0, system.box.yz )), box_unit) 
            yhi_bound = uc.get_in_units(system.box.yhi, box_unit) + uc.get_in_units(max(( 0.0, system.box.yz )), box_unit)
            zlo_bound = uc.get_in_units(system.box.zlo, box_unit)
            zhi_bound = uc.get_in_units(system.box.zhi, box_unit)
            
            f.write('%f %f %f\n' % ( xlo_bound, xhi_bound, uc.get_in_units(system.box.xy, box_unit) )) 
            f.write('%f %f %f\n' % ( ylo_bound, yhi_bound, uc.get_in_units(system.box.xz, box_unit) ))
            f.write('%f %f %f\n' % ( zlo_bound, zhi_bound, uc.get_in_units(system.box.yz, box_unit) ))

        #write atomic header info and prepare outarray for writing
        header = 'ITEM: ATOMS id'
        print_string = '%i'
        outarray = np.empty((system.natoms, len(prop_info['LAMMPS-attribute'])))
        start = 0
        for attr, a_keys in prop_info['LAMMPS-attribute'].iteritems():
            
            #get first prop relation for attr
            relation = a_keys.list('relation')[0]
            prop =   relation.get('prop')
            index = (Ellipsis, ) + tuple(relation.list('index'))
            unit =  relation.get('unit', None)
            if unit == 'scaled':
                unit = None
                scale = True
            else:
                scale = False            
            
            #pass values to outarray
            outarray[:,start] = uc.get_in_units(system.atoms_prop(key=prop, scale=scale), unit)[index].reshape((system.natoms))
            start += 1
            
            #prepare header and print_string
            header += ' %s' % attr
            if am.tools.is_dtype_int(system.atoms.dtype[prop]):
                print_string += ' %i'
            else:
                print_string += ' ' + xf
            
        f.write(header + '\n')
        print_string += '\n'
        
        #iterate over all atoms
        for i in xrange(system.natoms):
            vals = (i+1, ) + tuple(outarray[i])
            f.write(print_string % vals)
        

def __prop_info_default_load(name_list):
    """Construct default conversion model for atom_dump.load()"""
    name_list = list(name_list)
    prop_info = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate'] = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate']['box_prop'] = DataModelDict([('unit', None)])
    prop_info['LAMMPS-dump-atoms_prop-relate']['atoms_prop'] = props = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate']['LAMMPS-attribute'] = attrs = DataModelDict()
    
    props['atype'] = DataModelDict()
    props['atype']['dtype'] = 'int'  
    
    attrs['type'] = DataModelDict([('relation', DataModelDict([('prop', 'atype')]) )])
    
    props['pos'] = DataModelDict()
    props['pos']['dtype'] = 'float'
    props['pos']['shape'] = 3
    
    #Add x, y, z
    if 'x' in name_list and 'y' in name_list and 'z' in name_list:    
        attrs['x'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                ('index', 0)]) )])
        attrs['y'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                ('index', 1)]) )])
        attrs['z'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                ('index', 2)]) )])
        temp = name_list.pop(name_list.index('x'))
        temp = name_list.pop(name_list.index('y'))
        temp = name_list.pop(name_list.index('z'))
    
    #Add xs, ys, zs
    elif 'xs' in name_list and 'ys' in name_list and 'zs' in name_list:
        attrs['xs'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 0),
                                                                 ('unit', 'scaled')]) )])
        attrs['ys'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 1),
                                                                 ('unit', 'scaled')]) )])
        attrs['zs'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 2),
                                                                 ('unit', 'scaled')]) )])
        temp = name_list.pop(name_list.index('xs'))
        temp = name_list.pop(name_list.index('ys'))
        temp = name_list.pop(name_list.index('zs'))
    
    #Add xu, yu, zu    
    elif 'xu' in name_list and 'yu' in name_list and 'zu' in name_list:
        attrs['xu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 0)]) )])
        attrs['yu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 1)]) )])
        attrs['zu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                 ('index', 2)]) )])
        temp = name_list.pop(name_list.index('xu'))
        temp = name_list.pop(name_list.index('yu'))
        temp = name_list.pop(name_list.index('zu'))
        
    #Add xsu, ysu, zsu
    elif 'xsu' in name_list and 'ysu' in name_list and 'zsu' in name_list:
        attrs['xsu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                  ('index', 0),
                                                                  ('unit', 'scaled')]) )])
        attrs['ysu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                  ('index', 1),
                                                                  ('unit', 'scaled')]) )])
        attrs['zsu'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                                  ('index', 2),
                                                                  ('unit', 'scaled')]) )])
        temp = name_list.pop(name_list.index('xsu'))
        temp = name_list.pop(name_list.index('ysu'))
        temp = name_list.pop(name_list.index('zsu'))
        
    else:
        raise ValueError('Full 3D set of positions not in file')    
    
    for name in name_list:
        if name != 'id' and name != 'type':
            props[name] = DataModelDict()            
            attrs[name] = DataModelDict([('relation', DataModelDict([('prop', name)]) )])
    return prop_info        
            
def __prop_info_default_dump(system):
    """Construct default conversion model for atom_dump.load()"""

    #iterates over value to create unique parameter name
    def __itername(string, value):
        for i in xrange(len(value)):
            if hasattr(value[i], '__iter__'):
                down = __itername(string + '[%i]'%(i+1), value[i])
                for d in down:
                    yield d
            else:
                yield string + '[%i]'%(i+1)
    
    #iterates over value to return corresponding index values
    def __iterindex(tup, value):
        for i in xrange(len(value)):
            if hasattr(value[i], '__iter__'):
                down = __iterindex(tup + (i,), value[i])
                for d in down:
                    yield d
            else:
                yield tup + (i,)    
    
    prop_info = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate'] = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate']['box_prop'] = DataModelDict([('unit', None)])
    prop_info['LAMMPS-dump-atoms_prop-relate']['atoms_prop'] = props = DataModelDict()
    prop_info['LAMMPS-dump-atoms_prop-relate']['LAMMPS-attribute'] = attrs = DataModelDict() 
    
    props['atype'] = DataModelDict([('dtype', 'int')])
    
    attrs['type'] = DataModelDict([('relation', DataModelDict([('prop', 'atype')]) )])
    
    props['pos'] = DataModelDict([('dtype', 'float'),
                                  ('shape', 3)])
    
    attrs['x'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                            ('index', 0)]) )])
    attrs['y'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                            ('index', 1)]) )])
    attrs['z'] = DataModelDict([('relation', DataModelDict([('prop', 'pos'), 
                                                            ('index', 2)]) )])
    
    for prop in system.atoms_prop():
        if prop != 'atype' and prop != 'pos':
            p = system.atoms_prop(a_id=0, key=prop)
            
            #for scalars, no shape is defined and no index is needed to relate to
            if len(p.shape) == 0:
                props[prop] = DataModelDict()
                attrs[prop] = DataModelDict([('relation', DataModelDict([('prop', prop)]) )])
            
            #for vectors, reduce shape and index lists to single ints
            elif len(p.shape) == 1:
                props[prop] = DataModelDict([('shape', p.shape[0])])
                names = [val for val in __itername(prop, p)]
                indexes = [val for val in __iterindex((), p)]
                for name, index in zip(names, indexes):
                    attrs[name] = DataModelDict([('relation', DataModelDict([('prop', prop), 
                                                                             ('index', index[0])]) )])
            
            #for all others, include shape and indexes as lists
            else:
                props[prop] = DataModelDict([('shape', p.shape)])
                names = [val for val in __itername(prop, p)]
                indexes = [val for val in __iterindex((), p)]
                for name, index in zip(names, indexes):
                    attrs[name] = DataModelDict([('relation', DataModelDict([('prop', prop), 
                                                                             ('index', index)]) )])
    return prop_info        