# coding: utf-8

# atomman imports
from .. import load_system_model
from ...library import Database
from ...tools import screen_input

import pandas as pd

def load(key=None, method='dynamic', standing='good',
                   family=None, parent_key=None, potential=None,
                   potential_LAMMPS_id=None, potential_LAMMPS_key=None, potential_id=None, potential_key=None,
                   symbols=None, natoms=None, natypes=None, keyword=None,
                   verbose=False, database=None, localpath=None, local=False, remote=True):

    """
    Loads a potential-dependent relaxed crystal record from the library.  If
    multiple matches are found based on inputs a selection menu will appear.

    Parameters
    ----------
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    method : str or list or None, optional
        The relaxation method used.  Allowed values are dynamic, static and box.
        Default value is dynamic (the most rigorous relaxation method).  All
        will be loaded if set to None.
    standing : str or list or None, optional
        "good" records are the unique crystals found with the most rigorous
        relaxation, and with known prototypes over DFT structures.  "bad" are
        records filtered out, usually for being duplicates.  Default value is
        "good".  All will be loaded if set to None.
    family : str or atomman.library.CrystalPrototype or list, optional
        The crystal family associated with the relaxed crystal - either crystal
        prototype name or MP/OQMD database entry name.
    parent_key : str or list, optional
        The UUID4 key(s) assigned to the calculation that the record is based
        on.
    potential : atomman.lammps.Potential or list, optional
        A loaded LAMMPS potential object to limit the search by.
    potential_LAMMPS_id : str or list, optional
        The id for a LAMMPS implemented potential to limit the search by.
    potential_LAMMPS_key : str or list, optional
        The UUID4 for a LAMMPS implemented potential to limit the search by.
    potential_id : str or list, optional
        The id for a potential to limit the search by.
    potential_key : str or list, optional
        The UUID4 for a potential to limit the search by.
    symbols : str or list, optional
        Element symbols in the crystal to limit the search by.
    natypes : int or list, optional
        The number of unique element model symbols in the crystal to limit
        the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
    keyword : str, optional
        If given, will limit the search to all records that contain the keyword
        substring.  Cannot be combined with any of the above parameters.
    database : atomman.library.Database, optional
        A pre-defined Database object to use.  If not given, will initialize
        a new Database object.  Passing in a database can save time if multiple
        calls are made for the same record type. 
    localpath : str, optional
        The local library path to use when initializing a new Database.  IF not
        given, will use the default localpath.  Ignored if database is given. 
    local : bool, optional
        Indicates if the Database object is to look for local records.  Default
        is True.  Ignored if database is given.
    remote : bool, optional
        Indicates if the Database object is to look for remote records.  Default
        is True.  Ignored if database is given.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    system : atomman.System
        The system object generated from the relaxed crystal.
    """
    # Create Database object and load if needed
    if database is None:
        if local is True:
            database = Database(load='relaxed_crystal', localpath=localpath,
                                local=local, remote=remote, verbose=verbose)
        else:
            database = Database(local=local, remote=remote, verbose=verbose)
    
    # Get crystal prototype record
    crystal = database.get_relaxed_crystal(key=key, method=method, standing=standing,
                   family=family, parent_key=parent_key, potential=potential,
                   potential_LAMMPS_id=potential_LAMMPS_id, potential_LAMMPS_key=potential_LAMMPS_key,
                   potential_id=potential_id, potential_key=potential_key,
                   symbols=symbols, natoms=natoms, natypes=natypes, keyword=keyword,
                   verbose=verbose)
    
    # Retrieve unit cell information and set symbols
    return crystal.ucell