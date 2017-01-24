import atomman as am

import numpy as np

def nye_tensor_p(system, neighbor_list=None, neighbor_list_cutoff=None):
    """This is a utility function for generating the p vector list for more complex atomic systems"""

    #neighbor list setup
    if neighbor_list is not None:
        assert neighbor_list_cutoff is None, 'neighbor_list and neighbor_list_cutoff cannot both be given'
    elif neighbor_list_cutoff is not None:
        neighbor_list = am.nlist(system, neighbor_list_cutoff)
    elif 'nlist' in system.prop:
        neighbor_list = system.prop['nlist']    
    
    p = []
    for i in xrange(system.natoms):
        p.append(system.dvect(i, neighbor_list[i][1:neighbor_list[i][0]+1]))
    return np.asarray(p)