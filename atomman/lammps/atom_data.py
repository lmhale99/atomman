#External library imports
import numpy as np

#Internal library imports
import style
import atomman as am
import atomman.unitconvert as uc

def load(fname, pbc=(True, True, True), atom_style='atomic', units='metal'):
    """
    Read a LAMMPS-style atom data file and return a System.
    
    Argument:
    fname = name (and location) of file to read data from.
    
    Keyword Arguments:
    pbc -- list or tuple of three boolean values indicating which System directions are periodic. Default is (True, True, True).
    atom_style -- LAMMPS atom_style option associated with the data file.  Default is 'atomic'.
    units -- LAMMPS units option associated with the data file. Default is 'metal'.
    
    When the file is read in, the units of all property values are automatically converted to atomman's set working units.
    """
    
    units_dict = style.unit(units)

    readtime = False
    count = 0
    xy = 0.0
    xz = 0.0
    yz = 0.0
    system = None
    with open(fname, 'r') as fp:
        #loop over all lines in fp
        for line in fp:
            terms = line.split()
            if len(terms)>0:
                
                #read atomic information if time to do so
                if readtime == True:
                    a_id = int(terms[0]) - 1
                    prop_vals[a_id] = terms[1:]
                    count += 1
                    
                    #save values to system once all atoms read in
                    if count == natoms:
                        readtime = False
                        count = 0
                        start = 0
                        
                        #iterate over all atom_style properties
                        for name, v in props.iteritems():
                            if name != 'a_id':
                                size, dim, dtype = v
                                value = np.asarray(prop_vals[:, start:start+size], dtype=dtype)
                                start += size
                                
                                #set units according to LAMMPS units style
                                unit = units_dict.get(dim, None)
                                system.atoms_prop(key=name, value=uc.set_in_units(value, unit))
                
                #read number of atoms 
                elif len(terms) == 2 and terms[1] == 'atoms':
                    natoms = int(terms[0])
                
                #read number of atom types
                elif len(terms) == 3 and terms[1] == 'atom' and terms[2] == 'types': 
                    natypes = int(terms[0])
                
                #read boundary info
                elif len(terms) == 4 and terms[2] == 'xlo' and terms[3] == 'xhi':
                    xlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 4 and terms[2] == 'ylo' and terms[3] == 'yhi':
                    ylo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    yhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 4 and terms[2] == 'zlo' and terms[3] == 'zhi':
                    zlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    zhi = uc.set_in_units(float(terms[1]), units_dict['length'])
                elif len(terms) == 6 and terms[3] == 'xy' and terms[4] == 'xz' and terms[5] == 'yz':
                    xy = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xz = uc.set_in_units(float(terms[1]), units_dict['length'])  
                    yz = uc.set_in_units(float(terms[2]), units_dict['length'])
                
                #Flag when reached data and setup for reading
                elif len(terms) == 1 and terms[0] in ('Atoms', 'Velocities'):
                    #create system if not already
                    if system is None:
                        box = am.Box(xlo=xlo, xhi=xhi, 
                                    ylo=ylo, yhi=yhi, 
                                    zlo=zlo, zhi=zhi, 
                                    xy=xy, xz=xz, yz=yz)
                                  
                        system = am.System(box=box, atoms=am.Atoms(natoms=natoms), pbc = pbc)    
                    
                    if terms[0] == 'Atoms':
                        props = style.atom(atom_style)
                    else:
                        props = style.velocity(atom_style)
                    
                    nvals = 0
                    for name, v in props.iteritems():
                        nvals += v[0]
                    prop_vals = np.empty((natoms, nvals-1), dtype=float)
                    readtime = True
                     
    assert system.natypes == natypes, 'Number of atom types does not match!'
    return system

def dump(fname, system, units='metal', atom_style='atomic'):
    """
    Write a LAMMPS-style atom data file from a System.
    
    Argument:
    fname -- name (and location) of file to save data to.
    system -- System to write to the atom data file.
    
    Keyword Arguments:
    atom_style -- LAMMPS atom_style option associated with the data file.  Default is 'atomic'.
    units -- LAMMPS units option associated with the data file. Default is 'metal'.
    
    When the file is written, the units of all property values are automatically converted from atomman's working units to the LAMMPS style units.
    """
    #wrap atoms because LAMMPS hates atoms out of bounds in atom data files
    system.wrap()
    
    #get unit information according to the units style
    units_dict = style.unit(units)
    length = units_dict['length']
    
    #open file
    with open(fname,'w') as f:
    
        #header info
        f.write('\n%i atoms\n' % system.natoms)
        f.write('%i atom types\n' % system.natypes)
        f.write('%f %f xlo xhi\n' % (uc.get_in_units(system.box.xlo, length), uc.get_in_units(system.box.xhi, length)))
        f.write('%f %f ylo yhi\n' % (uc.get_in_units(system.box.ylo, length), uc.get_in_units(system.box.yhi, length)))
        f.write('%f %f zlo zhi\n' % (uc.get_in_units(system.box.zlo, length), uc.get_in_units(system.box.zhi, length)))
        
        #add tilts if tilt values not equal to zero
        if system.box.xy == 0.0 and system.box.xz == 0.0 and system.box.yz == 0.0:
            pass
        else:
            f.write('%f %f %f xy xz yz\n' % (uc.get_in_units(system.box.xy, length), uc.get_in_units(system.box.xz, length), uc.get_in_units(system.box.yz, length)))
        
        #Write atom info
        f.write('\nAtoms\n\n')
        props = style.atom(atom_style)
        
        #Count how many terms are being printed
        all_size = 0
        for v in props.itervalues():
            all_size += v[0]
        
        outarray = np.empty((system.natoms, all_size))
        start = 0
        print_string = ''
        
        #iterate over all properties and set values for printing
        for name, v in props.iteritems():
            size, dim, dtype = v
            try:
                unit = units_dict[dim]
            except:
                unit = None
            
            if dtype == int:
                for i in xrange(size):
                    print_string += ' %i'
            else:
                for i in xrange(size):  
                    print_string += ' %.13e'
            
            if name == 'a_id':
                outarray[:, start:start+size] = np.arange(1, system.natoms+1).reshape((system.natoms, size))
            else:
                outarray[:, start:start+size] = uc.get_in_units(system.atoms_prop(key=name), unit).reshape((system.natoms, size))
            start += size
        print_string = print_string.strip() + '\n'
        
        #iterate over all atoms
        for i in xrange(system.natoms):
            f.write(print_string % tuple(outarray[i]))
        
            
        #Test for velocity info
        v_test = system.atoms_prop(a_id=0, key='velocity')
        if v_test is not None:
            
            #Write velocity info
            f.write('\nVelocities\n\n')
            props = style.velocity(atom_style) 
        
            #Count how many terms are being printed
            all_size = 0
            for v in props.itervalues():
                all_size += v[0]
            
            outarray = np.empty((system.natoms, all_size))
            start = 0
            print_string = ''
            
            #iterate over all properties and set values for printing
            for name, v in props.iteritems():
                size, dim, dtype = v
                try:
                    unit = units_dict[dim]
                except:
                    unit = None

                if dtype == int:
                    for i in xrange(size):
                        print_string += ' %i'
                else:
                    for i in xrange(size):  
                        print_string += ' %.13e'
                
                if name == 'a_id':
                    outarray[:, start:start+size] = np.arange(1, system.natoms+1).reshape((system.natoms, size))
                else:
                    outarray[:, start:start+size] = uc.get_in_units(system.atoms_prop(key=name), unit).reshape((system.natoms, size))
                start += size
            print_string = print_string.strip() + '\n'
            
            #iterate over all atoms
            for i in xrange(system.natoms):
                f.write(print_string % tuple(outarray[i]))
    
    #return appropriate unts, atom_style, boundary, and read_data LAMMPS commands    
    boundary = ''
    for i in xrange(3):
        if system.pbc[i]:
            boundary += 'p '
        else:
            boundary += 'm '    
    
    newline = '\n'
    script = newline.join(['#Script and atom data file prepared by AtomMan package',
                           '',
                           'units ' + units,
                           'atom_style ' + atom_style,
                           ''
                           'boundary ' + boundary,
                           'read_data ' + fname])
    return script 