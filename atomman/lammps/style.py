from collections import OrderedDict
import numpy as np
    
def atom(atom_style='atomic'):
    """
    Returns the atom information associated with a given LAMMPS atom_style.
    
    Returns a dictionary where the keys are atom property names and values are tuples consisting of:
    (size, physical-quantity, data-type). 
    """
    int = np.dtype('int32')
    float = np.dtype('float64')
    
    params = OrderedDict()

    if atom_style == 'angle':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)        #integer molecule-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'atomic':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'body':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['b_flag'] =          (1, None, int)       #bodyflag
        params['mass'] =            (1, 'mass', float)      #particle mass
        params['pos'] =             (3, 'length', float)    #vector position
                
    elif atom_style == 'bond':
        params['a_id'] =            (1, None, int)       #integer atom-id
        params['m_id'] =            (1, None, int)       #integer molecule-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'charge':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['q'] =               (1, 'charge', float)    #atom charge
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'dipole':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['q'] =               (1, 'charge', float)    #atom charge
        params['pos'] =             (3, 'length', float)    #vector position
        params['mu'] =              (3, 'dipole', float)    #dipole moment
        
    elif atom_style == 'electron':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['q'] =               (1, 'charge', float)    #atom charge
        params['spin'] =            (1, None, int)       #electron spin
        params['e_radius'] =        (1, 'length', float)    #electron radius
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'ellipsoid':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['e_flag'] =          (1, None, int)       #ellipsoidflag
        params['density'] =         (1, 'density', float)   #particle density
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'full':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)        #integer molecule-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['q'] =               (1, 'charge', float)    #atom charge
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'line':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)        #integer molecule-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['l_flag'] =          (1, None, int)       #lineflag
        params['density'] =         (1, 'density', float)   #particle density
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'meso':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['rho'] =             (1, None, float)        #SPH density (not unit controlled?)
        params['e'] =               (1, None, float)        #SPH energy (not unit controlled?)
        params['cv'] =              (1, None, float)        #SPH heat capacity (not unit controlled?)
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'molecular':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)        #integer molecule-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'peri':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['volume'] =          (1, 'volume', float)    #particle volume
        params['density'] =         (1, 'density', float)   #particle density
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'smd':
        params['a_id'] =            (1, None, int)       #integer atom-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['m_id'] =            (1, None, int)       #integer molecule-id
        params['volume'] =          (1, 'volume', float)    #particle volume
        params['mass'] =            (1, 'mass', float)      #particle mass
        params['k_radius'] =        (1, 'length', float)    #kernal_radius
        params['c_radius'] =        (1, 'length', float)    #contact_radius
        params['pos'] =             (3, 'length', float)    #vector position
    
    elif atom_style == 'sphere':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)       #integer atom type
        params['diameter'] =        (1, 'length', float)   #particle diameter
        params['density'] =         (1, 'density', float)   #particle density
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'template':        
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)       #integer molecule-id
        params['m_temp'] =          (1, None, int)        #molecule template index
        params['a_temp'] =          (1, None, int)        #atom template index
        params['atype'] =           (1, None, int)        #integer atom type
        params['pos'] =             (3, 'length', float)    #vector position

    elif atom_style == 'tri':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['m_id'] =            (1, None, int)        #integer molecule-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['t_flag'] =          (1, None, int)       #lineflag
        params['density'] =         (1, 'density', float)   #particle density
        params['pos'] =             (3, 'length', float)    #vector position
        
    elif atom_style == 'wavepacket':
        params['a_id'] =            (1, None, int)        #integer atom-id
        params['atype'] =           (1, None, int)        #integer atom type
        params['q'] =               (1, 'charge', float)    #atom charge
        params['spin'] =            (1, None, int)        #electron spin
        params['e_radius'] =        (1, 'length', float)    #electron radius
        params['e_id'] =            (1, None, int)        #integer electron-id (etag)
        params['cs_re'] =           (1, None, float)        #real part of wavepacket coefficients
        params['cs_im'] =           (1, None, float)        #imaginary part of wavepacket coefficients
        params['pos'] =             (3, 'length', float)    #vector position      
    
    elif atom_style[:6] == 'hybrid':
        substyles = atom_style.split()
        params = atom('atomic')
        for substyle in substyles[1:]:
            subparams = atom(substyle)
            for k, v in subparams.iteritems():
                if k not in params:
                    params[k] = v
    else:
        raise ValueError('atom_style ' + atom_style + ' not supported')
    
    return params

def velocity(atom_style='atomic'):
    """
    Returns the velocity information associated with a given LAMMPS atom_style.
    
    Returns a dictionary where the keys are velocity property names and values are tuples consisting of:
    (size, physical-quantity, data-type). 
    """
    
    params = OrderedDict()
    
    if atom_style in ('angle', 'atomic', 'body', 'bond', 'charge', 'dipole', 'full', 'line', 'meso',
                      'molecular', 'peri', 'smd', 'template', 'tri', 'wavepacket'):
        params['a_id'] =            (1, None, int)            #integer atom-id
        params['velocity'] =        (3, 'velocity', float)      #translational velocity
   
    elif atom_style == 'electron':
        params['a_id'] =            (1, None, int)            #integer atom-id
        params['velocity'] =        (3, 'velocity', float)      #translational velocity
        params['erval'] =           (1, 'velocity', float)      #electron radial velocity
        
    elif atom_style == 'ellipsoid':
        params['a_id'] =            (1, None, int)            #integer atom-id
        params['velocity'] =        (3, 'velocity', float)      #translational velocity
        params['ang-momentum'] =    (3, 'ang-mom', float)       #components of angular momentum
        
    elif atom_style == 'sphere':
        params['a_id'] =            (1, None, int)            #integer atom-id
        params['velocity'] =        (3, 'velocity', float)      #translational velocity
        params['ang-velocity'] =    (3, 'ang-vel', float)       #components of angular velocity
        
    elif atom_style[:6] == 'hybrid':
        substyles = atom_style.split()
        params = velocity()
        for substyle in substyles[1:]:
            subparams = atom_style_params(atom_style=substyle)
            for k, v in subparams.iteritems():
                if k not in params:
                    params[k] = v
    else:
        raise ValueError('atom_style ' + atom_style + ' not supported')
    
    return params
    
def unit(units='metal'):
    """
    Returns the units information associated with a LAMMPS units option.
    
    Returns a dictionary where:
    keys are physical quantity types
    values are associated units for the LAMMPS units option specified.
    """
    params = OrderedDict()
    
    if units == 'lj':
        params['mass'] =                None
        params['length'] =              None
        params['time'] =                None
        params['energy'] =              None
        params['velocity'] =            None
        params['force'] =               None
        params['torque'] =              None
        params['temperature'] =         None
        params['pressure'] =            None
        params['dynamic viscosity'] =   None
        params['charge'] =              None
        params['dipole'] =              None
        params['electric field'] =      None
        params['density'] =             None
                  
    elif units == 'real':
        params['mass'] =                'g/mol'
        params['length'] =              'angstrom' 
        params['time'] =                'fs'
        params['energy'] =              'kcal/mol'
        params['velocity'] =            'angstrom/fs'
        params['force'] =               'kcal/(mol*angstrom)'
        params['torque'] =              'kcal/mol'
        params['temperature'] =         'K'
        params['pressure'] =            'atm'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              'e'
        params['dipole'] =              'e*angstrom'
        params['electric field'] =      'V/angstrom'
        params['density'] =             'g/cm^3'

    elif units == 'metal':
        params['mass'] =                'g/mol'
        params['length'] =              'angstrom'
        params['time'] =                'ps'
        params['energy'] =              'eV'
        params['velocity'] =            'angstrom/ps'
        params['force'] =               'eV/angstrom'
        params['torque'] =              'eV'
        params['temperature'] =         'K'
        params['pressure'] =            'bar'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              'e'
        params['dipole'] =              'e*angstrom'
        params['electric field'] =      'V/angstrom'
        params['density'] =             'g/cm^3'

    elif units == 'si':
        params['mass'] =                'kg'
        params['length'] =              'm'
        params['time'] =                's'
        params['energy'] =              'J'
        params['velocity'] =            'm/s'
        params['force'] =               'N'
        params['torque'] =              'N*m'
        params['temperature'] =         'K'
        params['pressure'] =            'Pa'
        params['dynamic viscosity'] =   'Pa*s'
        params['charge'] =              'C'
        params['dipole'] =              'C*m'
        params['electric field'] =      'V/m'
        params['density'] =             'kg/m^3'

    elif units == 'cgs':
        params['mass'] =                'g'
        params['length'] =              'cm'
        params['time'] =                's'
        params['energy'] =              'erg'
        params['velocity'] =            'cm/s'
        params['force'] =               'dyn'
        params['torque'] =              'dyn*cm'
        params['temperature'] =         'K'
        params['pressure'] =            '0.1*Pa'
        params['dynamic viscosity'] =   'Pa*s/10'
        params['charge'] =              '10*c0*C'
        params['dipole'] =              '10*c0*C*cm'
        params['electric field'] =      'c0*uV/cm'
        params['density'] =             'g/cm^3'

    elif units == 'electron':
        params['mass'] =                'amu'
        params['length'] =              'aBohr'
        params['time'] =                'fs'
        params['energy'] =              '2*Ry'
        params['velocity'] =            '2*Ry*aBohr/hbar'
        params['force'] =               '2*Ry/aBohr'
        params['torque'] =              'dyn*cm'
        params['temperature'] =         'K'
        params['pressure'] =            'Pa'
        params['charge'] =              'e'
        params['dipole'] =              '1e-21/c0*C*m'
        params['electric field'] =      'V/cm'
                        
    elif units == 'micro':
        params['mass'] =                'pg'
        params['length'] =              'um'
        params['time'] =                'us'
        params['energy'] =              'pg*um^2/us^2'
        params['velocity'] =            'um/us'
        params['force'] =               'pg*um/us^2'
        params['torque'] =              'pg*um^2/us^2'
        params['temperature'] =         'K'
        params['pressure'] =            'pg/(um*us^2)'
        params['dynamic viscosity'] =   'pg/(um*us)'
        params['charge'] =              '1e-12*C'
        params['dipole'] =              '1e-12*C*um'
        params['electric field'] =      'V/um'
        params['density'] =             'pg/um^3'

    elif units == 'nano':
        params['mass'] =                '1e-18*g'
        params['length'] =              'nm'
        params['time'] =                'ns'
        params['energy'] =              '1e-18*g*nm^2/ns^2'
        params['velocity'] =            'nm/ns'
        params['force'] =               '1e-18*g*nm/ns^2'
        params['torque'] =              '1e-18*g*nm^2/ns^2'
        params['temperature'] =         'K'
        params['pressure'] =            '1e-18*g/(nm*ns^2)'
        params['dynamic viscosity'] =   '1e-18*g/(nm*ns)'
        params['charge'] =              'e'
        params['dipole'] =              'e*nm'
        params['electric field'] =      'V/nm'
        params['density'] =             '1e-18*g/nm^3'

    else:
        raise ValueError('units ' + units + ' not supported')
    params['ang-mom'] = params['length'] + '*' + params['velocity'] + '*' + params['mass']
    params['ang-vel'] = '1/' + params['time']
    
    return params

def timestep(units='metal'):
    if   units == 'lj': return 0.005
    elif units == 'real': return 1.0
    elif units == 'metal': return 0.001
    elif units == 'si': return 1.0e-8
    elif units == 'cgs': return 1.0e-8
    elif units == 'electron': return 0.001
    elif units == 'micro': return 2.0
    elif units == 'nano': return 0.00045
    else:
        raise ValueError('units ' + units + ' not supported')






































