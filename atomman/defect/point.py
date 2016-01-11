import numpy as np
import atomman as am

def point(system, ptd_type=None, atype=None, pos=None, ptd_id=None, db_vect=None, scale=False):
    #Returns a new system containing a point defect. 
    
    #Check that ptd_type is a valid option
    assert ptd_type in ('v', 'i', 's', 'db'),    'Invalid ptd_type. Options are: v, i, s, or db'
    
    #If vacancy
    if ptd_type == 'v':
        assert atype is None,                   'atype is meaningless for vacancy insertion'
        assert db_vect is None,                 'db_vect is meaningless for vacancy insertion'
        return vacancy(system, pos=pos, ptd_id=ptd_id, scale=scale)
        
    #If interstitial
    elif ptd_type == 'i':
        assert ptd_id is None,                   'ptd_id is meaningless for interstitial insertion'
        assert db_vect is None,                  'db_vect is meaningless for interstitial insertion' 
        return interstitial(system, atype=atype, pos=pos, scale=scale)
                
    #If substitutional
    elif ptd_type == 's':
        assert db_vect is None,                  'db_vect is meaningless for substitutional insertion'
        return substitutional(system, atype=atype, pos=pos, ptd_id=ptd_id, scale=scale)
            
    #if dumbbell
    elif ptd_type == 'db':
        return dumbbell(system, atype=atype, pos=pos, ptd_id=ptd_id, db_vect=db_vect, scale=scale)
    
def vacancy(system, pos=None, ptd_id=None, scale=False):
    #Adds a vacancy to the system by specifying a position or atom id
    
    pos_list = system.atoms_view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if scale:
            pos = system.unscale(pos)
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    #create new system and copy values over
    d_system = am.System(box=system.box(), pbc=system.pbc(), atoms=am.Atoms(natoms=system.natoms()-1))
    for prop in system.atoms_prop():
        view = system.atoms_view[prop]
        view = np.vstack(( view[:ptd_id], view[ptd_id+1:] ))
        d_system.atoms(prop, view, dtype=system.atoms_dtype[prop])
    
    #add property old_id with each atom's original id
    if d_system.atoms(0, 'old_id') is None:
        d_system.atoms('old_id', np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms()) )), dtype=int)
    
    return d_system
          
def interstitial(system, atype=None, pos=None, scale=False):
    #Adds an interstitial atom to the system with specified atom type and position
  
    pos_list = system.atoms_view['pos']
    
    if scale:
        pos = system.unscale(pos)
    
    #Use isclose and where to check that no atoms are already at pos
    ptd_id = np.where(np.isclose(pos_list, pos).all(axis=1))
    assert len(ptd_id) == 1 and len(ptd_id[0]) == 0,                                'atom already at pos'
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    
    #create new system and copy values over
    d_system = am.System(box=system.box(), pbc=system.pbc(), atoms=am.Atoms(natoms=system.natoms()+1))
    for prop in system.atoms_prop():
        view = system.atoms_view[prop]
        view = np.vstack(( view, np.zeros_like(view[0]) ))
        d_system.atoms_prop(prop, view, dtype=system.atoms_dtype[prop])
    d_system.atoms(d_system.natoms()-1, 'atype', atype)
    d_system.atoms(d_system.natoms()-1, 'pos', pos)
    
    #add property old_id with each atom's original id
    if d_system.atoms(0, 'old_id') is None:
        d_system.atoms('old_id', np.arange(d_system.natoms()), dtype=int)
    else:
        old_id = max(system.atoms_view['old_id'])[0] + 1
        d_system.atoms(d_system.natoms()-1, 'old_id', old_id)
    
    return d_system
    
def substitutional(system, atype=None, pos=None, ptd_id=None, scale=False):
    #Adds a substitutional atom to the system by changing the atom type of an atom indicated by ptd_id or pos
    
    pos_list = system.atoms_view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if scale:
            pos = system.unscale(pos)
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    assert system.atoms(ptd_id, 'atype') != atype,                                  'identified atom is already of the specified atype'
    
    #create new system and copy values over
    d_system = am.System(box=system.box(), pbc=system.pbc(), atoms=am.Atoms(natoms=system.natoms()))
    for prop in system.atoms_prop():
        view = system.atoms_view[prop]
        view = np.vstack(( view[:ptd_id], view[ptd_id+1:], view[ptd_id] ))
        d_system.atoms(prop, view, dtype=system.atoms_dtype[prop])
    d_system.atoms(d_system.natoms()-1, 'atype', atype)
    
    #add property old_id with each atom's original id
    if d_system.atoms(0, 'old_id') is None:
        d_system.atoms('old_id', np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms()), ptd_id )), dtype=int)
        
    return d_system
        
def dumbbell(system, atype=None, pos=None, ptd_id=None, db_vect=None, scale=False):
    #Adds a dumbbell interstitial to the system with specified atom type, position, and dumbbell vector
    
    pos_list = system.atoms_view['pos']
    
    #if pos is supplied, use isclose and where to identify the id of the atom at pos
    if pos is not None:
        if scale:
            pos = system.unscale(pos)
        assert ptd_id is None,                                                      'pos and ptd_id cannot both be supplied'
        ptd_id = np.where(np.isclose(pos_list, pos).all(axis=1))
        assert len(ptd_id) == 1 and len(ptd_id[0]) == 1,                            'Unique atom at pos not identified'
        ptd_id = long(ptd_id[0][0])
    
    #test that ptd_id is a valid entry
    try:
        pos = pos_list[ptd_id]
    except:
        raise TypeError('Invalid ptd_id')
    
    assert isinstance(atype, (int, long)) and atype > 0,                            'atype must be a positive integer'   
    
    #create new system and copy values over
    d_system = am.System(box=system.box(), pbc=system.pbc(), atoms=am.Atoms(natoms=system.natoms()+1))
    for prop in system.atoms_prop():
        view = system.atoms_view[prop]
        view = np.vstack(( view[:ptd_id], view[ptd_id+1:], view[ptd_id], np.zeros_like(view[0]) ))
        d_system.atoms(prop, view, dtype=system.atoms_dtype[prop])
        
    d_system.atoms(d_system.natoms()-1, 'atype', atype)
    d_system.atoms(d_system.natoms()-2, 'pos', pos-db_vect)
    d_system.atoms(d_system.natoms()-1, 'pos', pos+db_vect)
    
    #add property old_id with each atom's original id
    if d_system.atoms(0, 'old_id') is None:
        d_system.atoms('old_id', np.hstack(( np.arange(0, ptd_id), np.arange(ptd_id+1, system.natoms()), ptd_id, system.natoms())), dtype=int)
    else:
        old_id = max(system.atoms_view['old_id'])[0] + 1
        d_system.atoms(d_system.natoms()-1, 'old_id', old_id)
    
    return d_system