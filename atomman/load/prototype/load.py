# coding: utf-8
 
# Standard Python imports
from __future__ import annotations
from typing import Optional, Union

# atomman imports
from ... import System
from ...library import Database

def load(name: Union[str, list, None] = None,
         id: Union[str, list, None] = None,
         key: Union[str, list, None] = None,
         commonname: Union[str, list, None] = None,
         prototype: Union[str, list, None] = None,
         pearson: Union[str, list, None] = None,
         strukturbericht: Union[str, list, None] = None,
         sg_number: Union[int, list, None] = None,
         sg_hm: Union[str, list, None] = None,
         sg_schoenflies: Union[str, list, None] = None,
         a: Optional[float] = None,
         b: Optional[float] = None,
         c: Optional[float] = None,
         alpha: Optional[float] = None,
         beta: Optional[float] = None,
         gamma: Optional[float] = None,
         symbols: Optional[tuple] = None,
         database: Optional[Database] = None,
         local: Optional[bool] = True,
         remote: Optional[bool] = True,
         prompt: bool = True,
         refresh_cache: bool = False,
         verbose: bool = False) -> System:
    """
    Loads a crystal prototype record from the library. If multiple matches
    are found based on inputs a selection menu will appear.
    
    name : str or list, optional.
        Record names to search for.  These should be the same values as id.
    id : str or list, optional
        Prototype ID(s) to search for.  These are unique identifiers for each
        prototype based on comm.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    commonname : str or list, optional
        Common name(s) to limit the search by.
    prototype : str or list, optional
        Prototype identifying composition(s) to limit the search by.
    pearson : str or list, optional
        The Pearson symbol(s) to limit the search by.
    strukturbericht : str or list, optional
        The strukturbericht identifier(s) to limit the search by.
    sg_number : int or list, optional
        The space group number(s) to limit the search by.
    sg_hm : str or list, optional
        The space group Hermann-Maguin identifier(s) to limit the search by.
    sg_schoenflies : str or list, optional
        The space group Schoenflies identifier(s) to limit the search by.
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
    database : atomman.library.Database, optional
        A pre-defined Database object to use.  If not given, will initialize
        a new Database object.  Passing in a database can save time if multiple
        calls are made for the same record type. 
    local : bool, optional
        Indicates if the Database object is to look for local records.  Default
        is True.  Ignored if database is given.
    remote : bool, optional
        Indicates if the Database object is to look for remote records.  Default
        is True.  Ignored if database is given.
    prompt : bool, optional
        If prompt=True (default) then a screen input will ask for a selection
        if multiple matching potentials are found.  If prompt=False, then an
        error will be thrown if multiple matches are found.
    refresh_cache : bool, optional
        If the local database is of style "local", indicates if the metadata
        cache file is to be refreshed.  If False,
        metadata for new records will be added but the old record metadata
        fields will not be updated.  If True, then the metadata for all
        records will be regenerated, which is needed to update the metadata
        for modified records.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    system : atomman.System
        The system object generated from the crystal prototype.
    """

    # Create Database object and load if needed
    if database is None:
        database = Database.Database()
    
    # Get crystal prototype record
    prototype = database.get_crystal_prototype(local=local, remote=remote, name=name,
                            id=id, key=key, commonname=commonname, prototype=prototype,
                            pearson=pearson, strukturbericht=strukturbericht,
                            sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
                            prompt=prompt, refresh_cache=refresh_cache, verbose=verbose)
    
    # Retrieve unit cell information and set symbols
    ucell = prototype.ucell
    ucell.symbols = symbols
    
    # Allow lattice constants to be set based on crystal family
    family = ucell.box.identifyfamily()
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