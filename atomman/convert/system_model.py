# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from collections import OrderedDict
import os
import sys

# http://www.numpy.org/
import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.core.Atoms
import atomman.core.Box
import atomman.core.System
import atomman.unitconvert as uc
from atomman.compatibility import range, iteritems
import atomman.crystal

def load(model, key='atomic-system', index=0):
    """
    Read in a data model containing a crystal structure.
    
    Parameters
    ----------
    model : str, file-like object or DataModelDict
        The data model to read.
    key : str, optional
        The key identifying the root element for the system definition.
        Default value is 'atomic-system'.
    index : int, optional.
        If the full model has multiple key entries, the index specifies which
        to access.  Default value is 0 (first, or only entry).
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    symbols : list
        The list of atomic symbols corresponding to the system's atom types.
        Will be a list of None if symbol information is not in the data model.
    """
    
    # Pull system model out of data model using key and index
    a_sys = DM(model).finds(key)
    if len(a_sys) == 0:
        raise KeyError(key + ' not found in model')
    try:
        a_sys = a_sys[index]
    except:
        raise IndexError('Invalid index ' + str(index) + ' for model key ' + key)
    
    # Extract crystal system and box values
    c_system = list(a_sys['cell'].keys())[0]
    cell = a_sys['cell'][c_system]
    
    if c_system == 'cubic': 
        a = b = c = uc.value_unit(cell['a'])
        alpha = beta = gamma = 90.0
    
    elif c_system == 'hexagonal':
        a = b = uc.value_unit(cell['a'])
        c = uc.value_unit(cell['c'])
        alpha = beta = 90.0
        gamma = 120.0
    
    elif c_system == 'tetragonal':
        a = b = uc.value_unit(cell['a'])
        c = uc.value_unit(cell['c'])
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
        beta = cell['beta']
        gamma = cell['gamma']
    
    box = atomman.core.Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    
    # Count atypes and generate list of symbols if given
    atoms = []
    scale = None
    
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
        symbols = [None for i in range(max(all_atypes))]
        
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
            
            for atype, symbol in iteritems(sym_dict):
                symbols[atype-1] = symbol
    
    # Read per-atom properties
    natoms = len(atypes)
    prop = OrderedDict()
    prop['atype'] = atypes
    prop['pos'] = np.zeros((natoms, 3), dtype='float64')
    count = 0
    
    pos_units = []
    for atom in a_sys.iteraslist('atom'):
        
        # Read in pos for atom and unit info
        prop['pos'][count] = uc.value_unit(atom['position'])
        pos_units.append(atom['position'].get('unit', None))
        
        # Add other per-atom properties
        for property in atom.iteraslist('property'):
            if property['name'] not in prop:
                value = uc.value_unit(property)
                shape = (natoms, ) + value.shape
                prop[property['name']] = np.zeros(shape, dtype=value.dtype)
                
            prop[property['name']][count] = uc.value_unit(property)
        count += 1
    
    pos_unit = np.unique(pos_units)
    assert len(pos_unit) == 1, 'Mixed units for positions'
    if pos_unit[0] == 'scaled':
        scale=True
    else:
        scale=False
    
    atoms = atomman.core.Atoms(**prop)
    system = atomman.core.System(box=box, atoms=atoms, scale=scale)
    
    return system, symbols

def dump(system, **kwargs):
    """
    Return a DataModelDict 'cell' representation of the system.
    
    Parameters
    ----------
    system : atomman.System
        The system to generate the data model for.
    box_unit : str, optional
        Length unit to use for the box. Default value is 'angstrom'.
    symbols : list, optional
        list of atom-model symbols corresponding to the atom types.
    elements : list, optional
        list of element tags corresponding to the atom types.
    prop_units : dict, optional
        dictionary where the keys are the property keys to include, and
        the values are units to use. If not given, only the positions in
        scaled units are included.
    a_std : float, optional
        Standard deviation of a lattice constant to include if available.
    b_std : float, optional
        Standard deviation of b lattice constant to include if available.
    c_std : float, optional
        Standard deviation of c lattice constant to include if available.
    """
    
    # Set default values
    box_unit = kwargs.get('box_unit', 'angstrom')
    
    symbols = kwargs.get('symbols', [None for i in range(system.natypes)])
    if not isinstance(symbols, list):
        symbols = [symbols]
    assert len(symbols) == system.natypes, 'Number of symbols does not match number of atom types'
    
    elements = kwargs.get('elements', [None for i in range(system.natypes)])
    if not isinstance(elements, list):
        elements = [elements]
    assert len(elements) == system.natypes, 'Number of elements does not match number of atom types'
    
    prop_units = kwargs.get('prop_units', {})
    if 'pos' not in prop_units:
        prop_units['pos'] = 'scaled'
    
    # Extract system values
    a = system.box.a
    b = system.box.b
    c = system.box.c
    alpha = system.box.alpha
    beta =  system.box.beta
    gamma = system.box.gamma
    
    if 'a_std' in kwargs and 'b_std' in kwargs and 'c_std' in kwargs:
        errors = True
        a_std = kwargs['a_std']
        b_std = kwargs['b_std']
        c_std = kwargs['c_std']
    else:
        errors = False
        a_std = None
        b_std = None
        c_std = None
    
    model = DM()
    model['cell'] = cell = DM()
    
    # Test crystal family
    c_family = atomman.crystal.identifyfamily(system.box)
    if c_family is None:
        c_family = 'triclinic'
    cell[c_family] = DM()
    
    if c_family == 'cubic':
        a_ave = (a + b + c) / 3
        if errors is True:
            a_std_ave = (a_std + b_std + c_std) / 3
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
    
    elif c_family == 'tetragonal':
        a_ave = (a + b) / 2
        if errors is True:
            a_std_ave = (a_std + b_std) / 2
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        
    elif c_family == 'orthorhombic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
    
    elif c_family == 'hexagonal':
        a_ave = (a + b) / 2
        if errors is True:
            a_std_ave = (a_std + b_std) / 2
        else:
            a_std_ave = None
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        
    elif c_family == 'rhombohedral':
        a_ave = (a + b + c) / 3
        alpha_ave = (alpha + beta + gamma) / 3
        if errors is True:
            a_std_ave = (a_std + b_std + c_std) / 3
        else:
            a_std_ave = None
        
        cell[c_family]['a'] = uc.model(a_ave, box_unit, error=a_std_ave)
        cell[c_family]['alpha'] = alpha_ave
        
    elif c_family == 'monoclinic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        cell[c_family]['beta'] = beta
        
    elif c_family == 'triclinic':
        cell[c_family]['a'] = uc.model(a, box_unit, error=a_std)
        cell[c_family]['b'] = uc.model(b, box_unit, error=b_std)
        cell[c_family]['c'] = uc.model(c, box_unit, error=c_std)
        cell[c_family]['alpha'] = alpha
        cell[c_family]['beta'] = beta
        cell[c_family]['gamma'] = gamma
    else:
        raise ValueError('Unknown crystal family')
    
    atype = system.atoms.atype
    aindex = atype - 1
    
    for i in range(system.natoms):
        atom = DM()
        
        atom['component'] = int(atype[i])

        symbol = symbols[aindex[i]]
        if symbol is not None:
            atom['symbol'] = symbol
            
        element = elements[aindex[i]]
        if element is not None:
            atom['element'] = element
        
        atom['position'] = DM()
        if prop_units['pos'] == 'scaled':
            atom['position']['value'] = list(system.atoms_prop(key='pos', scale=True)[i])
        else:
            atom['position']['value'] = list(uc.get_in_units(system.atoms_prop(key='pos')[i], prop_units['pos']))
        atom['position']['unit'] = prop_units['pos']
            
        for key, unit in iteritems(prop_units):
            if key != 'pos' and key != 'atype':
                value = uc.get_in_units(system.atoms_prop(key=key)[i], unit)
                try:
                    value = list(value)
                except:
                    pass
                prop = DM()
                prop['name'] = key 
                prop['value'] = value
                prop['unit'] =  unit
                
                atom.append('property', prop)
                
        model.append('atom', atom)
        
    return DM([('atomic-system', model)])
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            ucell, symbols = crystal(f)
        print('symbols =', symbols)
        print(ucell)