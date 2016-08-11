import atomman as am

def load(style, input, **kwargs):
    """Load a System."""

    if style == 'system_model':
        key = kwargs.get('key', 'atomic-system')
        index = kwargs.get('index', 0)
        system, symbols = am.convert.system_model.load(input, key, index)
        
    elif style == 'cif':
        data_set = kwargs.get('data_set', None)
        system, symbols = am.convert.cif.load(input)
        
    elif style == 'atom_data':
        pbc = kwargs.get('pbc', (True, True, True))
        atom_style = kwargs.get('atom_style', 'atomic')
        units = kwargs.get('units', 'metal')
        system = am.lammps.atom_data.load(input, pbc, atom_style, units)
        symbols = [None for i in xrange(system.natypes)]
        
    elif style == 'atom_dump':
        prop_info = kwargs.get('prop_info', None)
        system = am.lammps.atom_dump.load(input, prop_info)
        symbols = [None for i in xrange(system.natypes)]
        
    elif style == 'ase_Atoms':
        system, symbols = am.convert.ase_Atoms.load(input)
        
    elif style == 'pymatgen_Structure':
        system, symbols = am.convert.pymatgen_Structure.load(input)
    
    elif style == 'poscar':
        system, symbols = am.convert.poscar.load(input)
    
    else:
        raise ValueError('Unsupported style')
    
    return system, symbols