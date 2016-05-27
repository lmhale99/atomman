import numpy as np
import atomman as am
import atomman.unitconvert as uc

def point(system, ptd_type=None, atype=None, pos=None, ptd_id=None, db_vect=None, scale=False, atol=None):
    """
    Returns a new System where a point defect has been inserted.
    
    Keyword Arguments:
    system -- base System that the defect is added to.    
    ptd_type -- key indicating which type of defect to add:
                'v' -- vacancy.
                'i' -- positional interstitial.
                's' -- substitutional.
                'db' -- dumbbell interstitial.
    atype -- atom type for defect atom ('i', 's', 'db' styles).
    pos -- position for adding the defect atom (all styles).
    ptd_id -- atom id where defect is added.  Alternative to using pos ('v', 's', 'db' styles).
    db_vect -- vector associated with the dumbbell interstitial ('db' style).
    scale -- indicates if pos and db_vect are absolute (False) or box-relative (True). Default is False.
    atol -- absolute tolerance for position-based searching. Default is 1e-3 angstroms.
    
    Adds atom property old_id if it doesn't already exist that tracks the original atom ids
    """ 
    
    #Check that ptd_type is a valid option
    assert ptd_type in ('v', 'i', 's', 'db'),    'Invalid ptd_type. Options are: v, i, s, or db'
    
    #If vacancy
    if ptd_type == 'v':
        assert atype is None,                   'atype is meaningless for vacancy insertion'
        assert db_vect is None,                 'db_vect is meaningless for vacancy insertion'
        return vacancy(system, pos=pos, ptd_id=ptd_id, scale=scale, atol=atol)
        
    #If interstitial
    elif ptd_type == 'i':
        assert ptd_id is None,                   'ptd_id is meaningless for interstitial insertion'
        assert db_vect is None,                  'db_vect is meaningless for interstitial insertion' 
        return interstitial(system, atype=atype, pos=pos, scale=scale, atol=atol)
                
    #If substitutional
    elif ptd_type == 's':
        assert db_vect is None,                  'db_vect is meaningless for substitutional insertion'
        return substitutional(system, atype=atype, pos=pos, ptd_id=ptd_id, scale=scale, atol=atol)
            
    #if dumbbell
    elif ptd_type == 'db':
        return dumbbell(system, atype=atype, pos=pos, ptd_id=ptd_id, db_vect=db_vect, scale=scale, atol=atol)
    
def vacancy(system, pos=None, ptd_id=None, scale=False, atol=None):
    """
    Returns a new System where a vacancy point defect has been inserted.
    
    Keyword Arguments:
    system -- base System that the defect is added to.
    pos -- position of the atom to be removed.
    ptd_id -- id of the atom to be removed.  Alternative to using pos.
    scale -- if pos is given, indicates if pos is absolute (False) or box-relative (True). Default is False.
    
    Adds atom property old_id if it doesn't already exist that tracks the original atom ids.
    """ 
    
    pos_list = system.atoms.view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if atol is None: 
            atol = uc.set_in_units(1e-3, 'angstrom')
        if scale:
            pos = system.unscale(pos)
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos, atol=atol).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    #create new system and copy values over
    d_system = am.System(box=system.box, pbc=system.pbc, atoms=am.Atoms(natoms=system.natoms-1))
    for prop in system.atoms_prop():
        view = system.atoms.view[prop]
        value = np.asarray(np.vstack(( view[:ptd_id], view[ptd_id+1:] )), dtype=system.atoms.dtype[prop])
        d_system.atoms_prop(key=prop, value=value)
    
    #add property old_id with each atom's original id
    if d_system.atoms_prop(key='old_id') is None:
        d_system.atoms_prop(key='old_id', value=np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms) )), dtype='int32')
    
    return d_system
          
def interstitial(system, atype=None, pos=None, scale=False, atol=None):
    """
    Returns a new System where a positional interstitial point defect has been inserted.
    
    Keyword Arguments:
    system -- base System that the defect is added to.
    atype -- atom type for the interstitial atom.
    pos -- position for adding the interstitial atom.
    scale -- if pos is given, indicates if pos is absolute (False) or box-relative (True). Default is False.
    
    Adds atom property old_id if it doesn't already exist that tracks the original atom ids.
    """  
  
    pos_list = system.atoms.view['pos']
    
    if atol is None: 
        atol = uc.set_in_units(1e-3, 'angstrom')
    if scale:
        pos = system.unscale(pos)

    #Use isclose and where to check that no atoms are already at pos
    ptd_id = np.where(np.isclose(pos_list, pos, atol=atol).all(axis=1))
    assert len(ptd_id) == 1 and len(ptd_id[0]) == 0,                                'atom already at pos'
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    
    #create new system and copy values over
    d_system = am.System(box=system.box, pbc=system.pbc, atoms=am.Atoms(natoms=system.natoms+1))
    for prop in system.atoms_prop():
        view = system.atoms.view[prop]
        value = np.asarray(np.vstack(( view, np.zeros_like(view[0]) )), dtype=system.atoms.dtype[prop])
        d_system.atoms_prop(key=prop, value=value)
    d_system.atoms_prop(a_id=d_system.natoms-1, key='atype', value=atype)
    d_system.atoms_prop(a_id=d_system.natoms-1, key='pos',   value=pos)
    
    #add property old_id with each atom's original id
    if d_system.atoms_prop(key='old_id') is None:
        d_system.atoms_prop(key='old_id', value=np.arange(d_system.natoms), dtype='int32')
    else:
        old_id = max(system.atoms_prop(key='old_id')) + 1
        d_system.atoms_prop(a_id=d_system.natoms-1, key='old_id', value=old_id)
    
    return d_system
    
def substitutional(system, atype=None, pos=None, ptd_id=None, scale=False, atol=None):
    """
    Returns a new System where a substitutional point defect has been inserted.
    
    Keyword Arguments:
    system -- base System that the defect is added to.
    atype -- atom type to change the indicated atom to.
    pos -- position of the atom to be changed.
    ptd_id -- id of the atom to be changed.  Alternative to using pos.
    scale -- if pos is given, indicates if pos is absolute (False) or box-relative (True). Default is False.
    
    Adds atom property old_id if it doesn't already exist that tracks the original atom ids.
    """ 
    
    pos_list = system.atoms.view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if atol is None: 
            atol = uc.set_in_units(1e-3, 'angstrom')
        if scale:
            pos = system.unscale(pos)        
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos, atol=atol).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    assert system.atoms_prop(a_id=ptd_id, key='atype') != atype,                    'identified atom is already of the specified atype'
    
    #create new system and copy values over
    d_system = am.System(box=system.box, pbc=system.pbc, atoms=am.Atoms(natoms=system.natoms))
    for prop in system.atoms_prop():
        view = system.atoms.view[prop]
        value = np.asarray(np.vstack(( view[:ptd_id], view[ptd_id+1:], view[ptd_id] )), dtype=system.atoms.dtype[prop])
        d_system.atoms_prop(key=prop, value=value)
    d_system.atoms_prop(a_id=d_system.natoms-1, key='atype', value=atype)
    
    #add property old_id with each atom's original id
    if d_system.atoms_prop(key='old_id') is None:
        d_system.atoms_prop(key='old_id', value=np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms), ptd_id )), dtype='int32')
        
    return d_system
        
def dumbbell(system, atype=None, pos=None, ptd_id=None, db_vect=None, scale=False, atol=None):
    """
    Returns a new System where a dumbbell interstitial point defect has been inserted.
    
    Keyword Arguments:
    system -- base System that the defect is added to.    
    atype -- atom type for the atom in the dumbbell pair being added to the system.
    pos -- position of the system atom where the dumbbell pair is added.
    ptd_id -- id of the system atom where the dumbbell pair is added.  Alternative to using pos.
    db_vect -- vector associated with the dumbbell interstitial.
    scale -- indicates if pos and db_vect are absolute (False) or box-relative (True). Default is False.
    
    Adds atom property old_id if it doesn't already exist that tracks the original atom ids.
    """ 
    
    pos_list = system.atoms.view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if atol is None: 
            atol = uc.set_in_units(1e-3, 'angstrom')
        if scale:
            pos = system.unscale(pos)
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos, atol=atol).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    
    #unscale db_vect if scale is True
    if scale:
        db_vect = system.unscale(db_vect)
    
    #create new system and copy values over
    d_system = am.System(box=system.box, pbc=system.pbc, atoms=am.Atoms(natoms=system.natoms+1))
    for prop in system.atoms_prop():
        view = system.atoms.view[prop]
        value = np.asarray(np.vstack(( view[:ptd_id], view[ptd_id+1:], view[ptd_id], np.zeros_like(view[0]) )), dtype=system.atoms.dtype[prop])
        d_system.atoms_prop(key=prop, value=value)
        
    d_system.atoms_prop(a_id=d_system.natoms-1, key='atype', value=atype)
    d_system.atoms_prop(a_id=d_system.natoms-2, key='pos',   value=pos-db_vect)
    d_system.atoms_prop(a_id=d_system.natoms-1, key='pos',   value=pos+db_vect)
    
    #add property old_id with each atom's original id
    if d_system.atoms_prop(a_id=0, key='old_id') is None:
        d_system.atoms_prop(key='old_id', value=np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms), ptd_id, system.natoms)), dtype='int32')
    else:
        old_id = max(system.atoms_prop(key='old_id')) + 1
        d_system.atoms_prop(a_id=d_system.natoms-1, key='old_id', value=old_id)
    
    return d_system