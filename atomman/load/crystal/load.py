# coding: utf-8

# Standard Python imports
from __future__ import annotations
from typing import Optional, Union

# https://github.com/usnistgov/potentials
from potentials.record import Record

# atomman imports
from ... import System
from ...library import Database
from ...tools import aslist

def load(name: Union[str, list, None] = None,
         key: Union[str, list, None] = None,
         method: Union[str, list, None] = 'dynamic',
         standing: Union[str, list, None] = 'good',
         family: Union[str, list, None] = None,
         parent_key: Union[str, list, None] = None, 
         potential: Optional[Record] = None,
         potential_LAMMPS_id: Union[str, list, None] = None,
         potential_LAMMPS_key: Union[str, list, None] = None,
         potential_id: Union[str, list, None] = None,
         potential_key: Union[str, list, None] = None,
         symbols: Union[str, list, None] = None,
         natoms: Union[int, list, None] = None,
         natypes: Union[int, list, None] = None,
         prompt: bool = True,
         refresh_cache: bool = False,
         verbose: bool = False,
         database: Optional[Database] = None,
         local: bool = True,
         remote: bool = True) -> System:
    """
    Loads a potential-dependent relaxed crystal record from the library.  If
    multiple matches are found based on inputs a selection menu will appear.

    Parameters
    ----------
    name : str or list
        The record name(s) to parse by.  For relaxed crystal records, the
        names should correspond to the key.
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
    prompt : bool
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
        The system object generated from the relaxed crystal.
    """
    # Create Database object and load if needed
    if database is None:
        database = Database.Database()
    
    if potential is not None:
        try:
            assert potential_LAMMPS_id is None
            assert potential_LAMMPS_key is None
            assert potential_id is None
            assert potential_key is None
        except:
            raise ValueError('potential cannot be given with the other potential parameters')
        potential_LAMMPS_key = []
        for pot in aslist(potential):
            potential_LAMMPS_key.append(pot.key)

    # Get crystal prototype record
    crystal = database.get_relaxed_crystal(local=local, remote=remote, name=name,
                           key=key, method=method, standing=standing, family=family,
                           parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
                           potential_LAMMPS_key=potential_LAMMPS_key,
                           potential_id=potential_id, potential_key=potential_key,
                           symbols=symbols, natoms=natoms, natypes=natypes,
                           prompt=prompt, refresh_cache=refresh_cache, verbose=verbose)
    
    # Retrieve unit cell information and set symbols
    return crystal.ucell