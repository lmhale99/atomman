import atomman as am

import numpy as np

def slip_vector(system_0, system_1, neighbor_list=None, neighbor_list_cutoff=None):
    """Compute the slip vectors for all atoms in system_1 relative to system_0."""

    assert system_0.natoms == system_1.natoms,  'systems have different number of atoms'

    #neighbor list setup
    if neighbor_list is not None:
        assert neighbor_list_cutoff is None, 'neighbor_list and neighbor_list_cutoff cannot both be given'
    elif neighbor_list_cutoff is not None:
        neighbor_list = am.nlist(system_0, neighbor_list_cutoff)
    elif 'nlist' in system_0.prop:
        neighbor_list = system_0.prop['nlist']
    
    #Calculate the slip vector
    slip = np.zeros((system_0.natoms, 3))
    
    for i in xrange(system_0.natoms):
        js = neighbor_list[i, 1:neighbor_list[i, 0]+1]
        slip[i] = -np.sum(system_1.dvect(i, js) - system_0.dvect(i, js), axis=0)
            
    return slip
    