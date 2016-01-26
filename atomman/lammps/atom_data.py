#External library imports
import numpy as np

#Internal library imports
import style
from atomman import Atoms, Box, System
from atomman.tools import unitconvert as uc

def load(fname, pbc=(True, True, True), atom_style='atomic', units='metal'):
    #Reads a LAMMPS-style atom data file and returns an atomman.System.
    
    units_dict = style.unit(units)

    with open(fname, 'r') as fin:
        readtime = False
        count = 0
        xy = None
        xz = None
        yz = None
        sys = None
        
        #loop over all lines in file
        for line in fin:
            terms = line.split()
            if len(terms)>0:
                
                #read atomic information if time to do so
                if readtime == True:
                    a_id = int(terms[0]) - 1
                    prop_vals[a_id] = terms[1:]
                    count += 1
                    
                    #save values to sys once all atoms read in
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
                                try:
                                    unit = units_dict[dim]
                                except:
                                    unit = None

                                sys.atoms(name, value, unit=unit)
                
                #read number of atoms 
                elif len(terms) == 2 and terms[1] == 'atoms':
                    natoms = int(terms[0])
                
                #read number of atom types
                elif len(terms) == 3 and terms[1] == 'atom' and terms[2] == 'types': 
                    natypes = int(terms[0])
                
                #read boundary info
                elif len(terms) == 4 and terms[2] == 'xlo' and terms[3] == 'xhi':
                    xlo = float(terms[0])
                    xhi = float(terms[1])
                elif len(terms) == 4 and terms[2] == 'ylo' and terms[3] == 'yhi':
                    ylo = float(terms[0])
                    yhi = float(terms[1])
                elif len(terms) == 4 and terms[2] == 'zlo' and terms[3] == 'zhi':
                    zlo = float(terms[0])
                    zhi = float(terms[1]) 
                elif len(terms) == 6 and terms[3] == 'xy' and terms[4] == 'xz' and terms[5] == 'yz':
                    xy = float(terms[0])
                    xz = float(terms[1])   
                    yz = float(terms[2])
                
                #Flag when reached data and setup for reading
                elif len(terms) == 1 and terms[0] in ('Atoms', 'Velocities'):
                    #create sys if not already
                    if sys is None:
                        box = Box(xlo=xlo, xhi=xhi, 
                                  ylo=ylo, yhi=yhi, 
                                  zlo=zlo, zhi=zhi, 
                                  xy=xy, xz=xz, yz=yz, 
                                  unit=units_dict['length'])
                                  
                        sys = System(box=box, 
                                     atoms=Atoms(natoms=natoms),
                                     pbc = pbc)    
                    
                    if terms[0] == 'Atoms':
                        props = style.atom(atom_style)
                    else:
                        props = style.velocity(atom_style)
                    
                    nvals = 0
                    for name, v in props.iteritems():
                        nvals += v[0]
                    prop_vals = np.empty((natoms, nvals-1), dtype=float)
                    readtime = True
                     
    assert sys.natypes() == natypes, 'Number of atom types does not match!'
    return sys 

def dump(fname, system, units='metal', atom_style='atomic'):
    #Writes a LAMMPS-style atom data file from a system. Returns LAMMPS commands for reading in the data file.

    #wrap atoms because LAMMPS hates atoms out of bounds in atom data files
    system.wrap()
    
    #get unit information according to the units style
    units_dict = style.unit(units)
    length = units_dict['length']
    
    #open file
    with open(fname,'w') as f:
    
        #header info
        f.write('\n%i atoms\n' % system.natoms())
        f.write('%i atom types\n' % system.natypes())
        f.write('%f %f xlo xhi\n' % (system.box('xlo', length), system.box('xhi', length)))
        f.write('%f %f ylo yhi\n' % (system.box('ylo', length), system.box('yhi', length)))
        f.write('%f %f zlo zhi\n' % (system.box('zlo', length), system.box('zhi', length)))
        
        #add tilts if tilt values not equal to zero
        if system.box('xy') == 0.0 and system.box('xz') == 0.0 and system.box('yz') == 0.0:
            pass
        else:
            f.write('%f %f %f xy xz yz\n' % (system.box('xy', length), system.box('xz', length), system.box('yz', length)))
        
        #Write atom info
        f.write('\nAtoms\n\n')
        props = style.atom(atom_style)
        
        #iterate over all atoms
        for i in xrange(system.natoms()):
            line = ''
            
            #iterate over all atom_style properties
            for name, v in props.iteritems():
                size, dim, dtype = v
                
                #use atom index + 1 for a_id
                if name == 'a_id':
                    value = i+1
                
                #get value in appropriate units
                else:
                    try:
                        unit = units_dict[dim]
                    except:
                        unit = None
                    value = system.atoms(i, name, unit=unit)                        

                #add property values to the line 
                if isinttype(value):
                    try:
                        for val in value.flatten():
                            line += '%i ' % val
                    except:
                        line += '%i ' % value
                else:
                    try:
                        for val in value.flatten():
                            line += '%.13e ' % val
                    except:
                        line += '%.13e ' % value
                
            f.write(line.strip()+'\n')
        
            
        #Test for velocity info
        v_test = system.atoms(0, 'velocity')
        if v_test is not None:
            
            #Write velocity info
            f.write('\nVelocities\n\n')
            props = style.velocity(atom_style) 
        
            #iterate over all atoms
            for i in xrange(system.natoms()):
                line = ''
                
                #iterate over all atom_style properties
                for name, v in props.iteritems():
                    size, dim, dtype = v
                    
                    #use atom index + 1 for a_id
                    if name == 'a_id':
                        value = i+1
                    
                    #get value in appropriate units
                    else:
                        try:
                            unit = units_dict[dim]
                        except:
                            unit = None
                        value = system.atoms(i, name, unit=unit)
                    
                    #add property values to the line 
                    if isinttype(value):
                        try:
                            for val in value.flatten():
                                line += '%i ' % val
                        except:
                            line += '%i ' % value
                    else:
                        try:
                            for val in value.flatten():
                                line += '%.13e ' % val
                        except:
                            line += '%.13e ' % value
                    
                f.write(line.strip()+'\n')
    
    #return appropriate unts, atom_style, boundary, and read_data LAMMPS commands    
    boundary = ''
    for i in xrange(3):
        if system.pbc(i):
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

def isinttype(value):
    if isinstance(value, (int, long)): 
        return True
    elif isinstance(value, np.ndarray):
        if value.dtype == int or value.dtype == long: 
            return True
        else:
            try:
                if issubclass(value.dtype.type, np.int):
                    return True
                else:
                    return False
            except:
                return False
    else:
        return False 
 