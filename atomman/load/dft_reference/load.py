# coding: utf-8

# Standard Python imports
from __future__ import annotations
from typing import Optional

# atomman imports
from ... import System
from ...library import Database

def load(id: str,
         api_key: Optional[str] = None,
         database: Optional[Database] = None,
         local: Optional[bool] = True,
         remote: Optional[bool] = True, 
         refresh_cache: bool = False,
         verbose: bool = False) -> System:
    """
    Loads a reference crystal structure from a DFT database.  Will first search
    the local database, then the remote NIST database to see if a matching
    entry is already included.  If a match is not found, then the structure
    will be obtained from the DFT database

    Parameters
    ----------
    id : str
        The reference crystal's unique id.  Combines a database tag "mp-" or
        "oqmd-" and the DFT database's entry id.
    api_key : str, optional
        The user's Materials Project API key or path to a file containing the
        key. Only needed for fetching structures from Materials Project and if
        the key is not set to the "MAPI_KEY" environment variable.
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
        The system object generated from the reference crystal.
    """
    # Create Database object and load if needed
    if database is None:
        database = Database.Database()
    
    # Fetch crystal from NIST or DFT database
    crystal = database.fetch_reference_crystal(id, api_key=api_key,
                                               local=local,
                                               remote=remote, 
                                               refresh_cache=refresh_cache,
                                               verbose=verbose)
    
    return crystal.ucell