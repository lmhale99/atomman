# coding: utf-8
 
# atomman imports
from .. import load_system_model
from ... import Library
from ...tools import identifyfamily, screen_input

def load(keyword=None, title=None, a=None, b=None, c=None,
         alpha=None, beta=None, gamma=None, symbols=None, local=True,
         remote=True):
    """
    Loads a crystal prototype record from the library. If multiple matches
    are found based on inputs a selection menu will appear.
    
    Parameters
    ----------
    keyword : str, optional
        Used to parse which prototype to load: a substring search is performed
        on the contents of all crystal prototype records and any containing
        this string are selected.
    title : str, optional
        The title (id) of the prototype record to load. Cannot be used with
        keyword.
    a : float, optional
        The a lattice parameter to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    b : float, optional
        The b lattice parameter to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    c : float, optional
        The c lattice parameter to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    alpha : float, optional
        The alpha lattice angle to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    beta : float, optional
        The beta lattice angle to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    gamma : gamma, optional
        The alpha lattice angle to scale the prototype by. Can only be given
        if it is a unique lattice parameter for the prototype's crystal family,
        and if all other unique lattice parameters are given.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
        Useful if the symbols for the model differ from the standard element
        tags or if the poscar file has no elemental information.
    local : bool, optional
        Indicates if the local version of the library should be searched.
        Default value is True.
    remote : bool, optional
        Indicates if the remote version of the library should be searched.
        Default value is True.
    
    Returns
    -------
    system : atomman.System
        The system object generated from the crystal prototype.
    """

    # Load library and fetch matching crystal_prototype records
    library = Library(remote=remote, local=local)
    records = library.potdb.get_records(template='crystal_prototype',
                                        keyword=keyword, title=title)

    # Check number of matches and select
    if len(records) == 1:
        record = records[0]
    elif len(records) > 1:
        print('Multiple matching prototypes found.')
        for i, record in enumerate(records):
            print(i+1, record['crystal-prototype']['id'])
        choice = int(screen_input('Select which one:'))
        record = records[choice-1]        
    else:
        raise ValueError('No matching prototypes found')
    
    ucell = load_system_model(record, symbols=symbols)
    
    family = identifyfamily(ucell.box)
    if family == 'cubic':
        if a is not None:
            if b is not None or c is not None or alpha is not None or beta is not None or gamma is not None:
                raise ValueError('Only a can be set for cubic prototypes')
            ucell.box_set(a=a, b=a, c=a, scale=True)
            
    elif family == 'hexagonal':
        if a is not None and c is not None:
            if b is not None or alpha is not None or beta is not None or gamma is not None:
                raise ValueError('Only a, c can be set for hexagonal prototypes')
            ucell.box_set(a=a, b=a, c=c, gamma=120, scale=True)
        elif a is not None or c is not None:
            raise ValueError('All or neither of a, c must be set for hexagonal prototypes')
            
    elif family == 'tetragonal':
        if a is not None and c is not None:
            if b is not None or alpha is not None or beta is not None or gamma is not None:
                raise ValueError('Only a, c can be set for tetragonal prototypes')
            ucell.box_set(a=a, b=a, c=c, scale=True)
        elif a is not None or c is not None:
            raise ValueError('All or neither of a, c must be set for tetragonal prototypes')
            
    elif family == 'rhombohedral':
        if a is not None and alpha is not None:
            if b is not None or c is not None or beta is not None or gamma is not None:
                raise ValueError('Only a, alpha can be set for rhombohedral prototypes')
            ucell.box_set(a=a, b=a, c=a, alpha=alpha, beta=alpha, gamma=alpha, scale=True)
        elif a is not None or alpha is not None:
            raise ValueError('All or neither of a, alpha must be set for rhombohedral prototypes')
            
    elif family == 'orthorhombic':
        if a is not None and b is not None and c is not None:
            if alpha is not None or beta is not None or gamma is not None:
                raise ValueError('Only a, b, c can be set for orthorhombic prototypes')
            ucell.box_set(a=a, b=b, c=c, scale=True)
        elif a is not None or b is not None or c is not None:
            raise ValueError('All or neither of a, b, c must be set for orthorhombic prototypes')
            
    elif family == 'monoclinic':
        if a is not None and b is not None and c is not None and beta is not None:
            if alpha is not None or gamma is not None:
                raise ValueError('Only a, b, c, beta can be set for monoclinic prototypes')
            ucell.box_set(a=a, b=b, c=c, beta=beta, scale=True)
        elif a is not None or b is not None or c is not None or beta is not None:
            raise ValueError('All or neither of a, b, c, beta must be set for monoclinic prototypes')
            
    else:
        if a is not None and b is not None and c is not None and alpha is not None and beta is not None and gamma is not None:
            ucell.box_set(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, scale=True)
        elif a is not None or b is not None or c is not None or alpha is not None or beta is not None or gamma is not None:
            raise ValueError('All or neither of a, b, c, alpha, beta, gamma must be set for triclinic prototypes')
            
    return ucell