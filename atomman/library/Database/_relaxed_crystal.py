# coding: utf-8

# Standard Python imports
from typing import Optional, Union, Tuple
from pathlib import Path

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

def get_relaxed_crystals(self,
                         name: Union[str, list, None] = None,
                         key: Union[str, list, None] = None,
                         method: Union[str, list, None] = None,
                         standing: Union[str, list, None] = None,
                         family: Union[str, list, None] = None,
                         parent_key: Union[str, list, None] = None, 
                         potential_LAMMPS_id: Union[str, list, None] = None,
                         potential_LAMMPS_key: Union[str, list, None] = None,
                         potential_id: Union[str, list, None] = None,
                         potential_key: Union[str, list, None] = None,
                         symbols: Union[str, list, None] = None,
                         natoms: Union[int, list, None] = None,
                         natypes: Union[int, list, None] = None,
                         local: Optional[bool] = None,
                         remote: Optional[bool] = None,
                         refresh_cache: bool = False,
                         return_df: bool = False,
                         verbose: bool = False
                         ) -> Union[np.ndarray, Tuple[np.ndarray, pd.DataFrame]]:
    """
    Gets all matching relaxed crystals from the database.
    
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
    standing : str or list or None, optional
        "good" records are the unique crystals found with the most rigorous
        relaxation, and with known prototypes over DFT structures.  "bad" are
        records filtered out, usually for being duplicates.
    family : str or atomman.library.CrystalPrototype or list, optional
        The crystal family associated with the relaxed crystal - either crystal
        prototype name or MP/OQMD database entry name.
    parent_key : str or list, optional
        The UUID4 key(s) assigned to the calculation that the record is based
        on.
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
        The number(s) of unique atom types to limit the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    refresh_cache : bool, optional
        If the local database is of style "local", indicates if the metadata
        cache file is to be refreshed.  If False,
        metadata for new records will be added but the old record metadata
        fields will not be updated.  If True, then the metadata for all
        records will be regenerated, which is needed to update the metadata
        for modified records.
    return_df : bool, optional
        If True, then the corresponding pandas.Dataframe of metadata
        will also be returned.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.

    Returns
    -------
    numpy.ndarray
        The matching record objects
    pandas.DataFrame
        A table of the associated record metadata, returned if return_df = True.
    """

    return self.get_records(
        style='relaxed_crystal', name=name,local=local, remote=remote, 
        refresh_cache=refresh_cache, return_df=return_df, verbose=verbose,
        key=key, method=method, standing=standing, family=family,
        parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
        potential_LAMMPS_key=potential_LAMMPS_key,
        potential_id=potential_id, potential_key=potential_key,
        symbols=symbols, natoms=natoms, natypes=natypes)

def promptfxn(df):
    """
    promptfxn for relaxed_crystal records
    """
    header = '#  family               symbols  alat    Ecoh    method  standing'
    print(header)

    js = df.sort_values('cohesive_energy').index
    for i, j in enumerate(js):
        crystal = df.loc[j]
        row =  f'{i+1:2} {crystal.family:20.20} '
        row += f'{"".join(crystal.symbols):8.8} '
        row += f'{crystal.a:7.4f} '
        row += f'{crystal.cohesive_energy:7.4f} '
        row += f'{crystal.method:7.7} '
        row += f'{crystal.standing:4.4}'
        print(row)
    
    i = int(input('Please select one:')) - 1
    if i < 0 or i >= len(js):
        raise ValueError('Invalid selection')

    return js[i]

def get_relaxed_crystal(self,
                        name: Union[str, list, None] = None,
                        key: Union[str, list, None] = None,
                        method: Union[str, list, None] = None,
                        standing: Union[str, list, None] = None,
                        family: Union[str, list, None] = None,
                        parent_key: Union[str, list, None] = None, 
                        potential_LAMMPS_id: Union[str, list, None] = None,
                        potential_LAMMPS_key: Union[str, list, None] = None,
                        potential_id: Union[str, list, None] = None,
                        potential_key: Union[str, list, None] = None,
                        symbols: Union[str, list, None] = None,
                        natoms: Union[int, list, None] = None,
                        natypes: Union[int, list, None] = None,
                        local: Optional[bool] = None,
                        remote: Optional[bool] = None,
                        prompt: bool = True,
                        refresh_cache: bool = False,
                        verbose: bool = False) -> Record:
    """
    Gets exactly one matching relaxed crystal from the database.
    
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
    standing : str or list or None, optional
        "good" records are the unique crystals found with the most rigorous
        relaxation, and with known prototypes over DFT structures.  "bad" are
        records filtered out, usually for being duplicates.
    family : str or atomman.library.CrystalPrototype or list, optional
        The crystal family associated with the relaxed crystal - either crystal
        prototype name or MP/OQMD database entry name.
    parent_key : str or list, optional
        The UUID4 key(s) assigned to the calculation that the record is based
        on.
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
    natypes : int or list, optional
        The number(s) of unique atom types to limit the search by.
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
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
    yabadaba.record.Record
        The matching record object
    """

    return self.get_record(
        style='relaxed_crystal', name=name, local=local, remote=remote, prompt=prompt,
        promptfxn=promptfxn, refresh_cache=refresh_cache, verbose=verbose,
        key=key, method=method, standing=standing, family=family,
        parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
        potential_LAMMPS_key=potential_LAMMPS_key,
        potential_id=potential_id, potential_key=potential_key,
        symbols=symbols, natoms=natoms, natypes=natypes)

def retrieve_relaxed_crystal(self,
                             name: Union[str, list, None] = None,
                             dest: Optional[Path] = None,
                             key: Union[str, list, None] = None,
                             method: Union[str, list, None] = None,
                             standing: Union[str, list, None] = None,
                             family: Union[str, list, None] = None,
                             parent_key: Union[str, list, None] = None, 
                             potential_LAMMPS_id: Union[str, list, None] = None,
                             potential_LAMMPS_key: Union[str, list, None] = None,
                             potential_id: Union[str, list, None] = None,
                             potential_key: Union[str, list, None] = None,
                             symbols: Union[str, list, None] = None,
                             natoms: Union[int, list, None] = None,
                             natypes: Union[int, list, None] = None,
                             local: Optional[bool] = None,
                             remote: Optional[bool] = None,
                             prompt: bool = True,
                             format: str = 'json',
                             indent: int = 4,
                             refresh_cache: bool = False,
                             verbose: bool = False):
    """
    Gets a single matching relaxed crystal from the database and saves it to a
    file based on the record's name.

    Parameters
    ----------
    name : str or list, optional
        The name(s) of records to limit the search by.
    dest : path, optional
        The parent directory where the record will be saved to.  If not given,
        will use the current working directory.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    method : str or list or None, optional
        The relaxation method used.  Allowed values are dynamic, static and box.
    standing : str or list or None, optional
        "good" records are the unique crystals found with the most rigorous
        relaxation, and with known prototypes over DFT structures.  "bad" are
        records filtered out, usually for being duplicates.
    family : str or atomman.library.CrystalPrototype or list, optional
        The crystal family associated with the relaxed crystal - either crystal
        prototype name or MP/OQMD database entry name.
    parent_key : str or list, optional
        The UUID4 key(s) assigned to the calculation that the record is based
        on.
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
    natypes : int or list, optional
        The number(s) of unique atom types to limit the search by.
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    prompt : bool, optional
        If prompt=True (default) then a screen input will ask for a selection
        if multiple matching potentials are found.  If prompt=False, then an
        error will be thrown if multiple matches are found.
    format : str, optional
        The file format to save the record in: 'json' or 'xml'.  Default
        is 'json'.
    indent : int, optional
        The number of space indentation spacings to use in the saved
        record for the different tiered levels.  Default is 4.  Giving None
        will create a compact record.
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
    
    Raises
    ------
    ValueError
        If local or remote is set to True when the corresponding database
        interaction has not been set.
    ValueError
        If multiple or no matching records are discovered.
    """
    self.retrieve_record(
        style='relaxed_crystal', name=name, dest=dest, local=local,
        remote=remote, prompt=prompt, promptfxn=promptfxn, format=format,
        indent=indent, refresh_cache=refresh_cache, verbose=verbose,
        key=key, method=method, standing=standing, family=family,
        parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
        potential_LAMMPS_key=potential_LAMMPS_key,
        potential_id=potential_id, potential_key=potential_key,
        symbols=symbols, natoms=natoms, natypes=natypes)

def download_relaxed_crystals(self, 
                              name: Union[str, list, None] = None,
                              key: Union[str, list, None] = None,
                              method: Union[str, list, None] = None,
                              standing: Union[str, list, None] = None,
                              family: Union[str, list, None] = None,
                              parent_key: Union[str, list, None] = None, 
                              potential_LAMMPS_id: Union[str, list, None] = None,
                              potential_LAMMPS_key: Union[str, list, None] = None,
                              potential_id: Union[str, list, None] = None,
                              potential_key: Union[str, list, None] = None,
                              symbols: Union[str, list, None] = None,
                              natoms: Union[int, list, None] = None,
                              natypes: Union[int, list, None] = None,
                              overwrite: bool = False,
                              return_records: bool = False,
                              verbose: bool = False) -> Optional[np.ndarray]:
    """
    Download citation records from the remote and save to localpath.
    
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
    standing : str or list or None, optional
        "good" records are the unique crystals found with the most rigorous
        relaxation, and with known prototypes over DFT structures.  "bad" are
        records filtered out, usually for being duplicates.
    family : str or atomman.library.CrystalPrototype or list, optional
        The crystal family associated with the relaxed crystal - either crystal
        prototype name or MP/OQMD database entry name.
    parent_key : str or list, optional
        The UUID4 key(s) assigned to the calculation that the record is based
        on.
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
    overwrite : bool, optional
        Flag indicating if any existing local records with names matching
        remote records are updated (True) or left unchanged (False).  Default
        value is False.
    return_records : bool, optional
        If True, the retrieved record objects are also returned.  Default
        value is False.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    numpy.ndarray
        The matching records - returned if return_records=True
    """

    return self.download_records(
        style='relaxed_crystal', name=name, overwrite=overwrite,
        return_records=return_records, verbose=verbose,
        key=key, method=method, standing=standing, family=family,
        parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
        potential_LAMMPS_key=potential_LAMMPS_key,
        potential_id=potential_id, potential_key=potential_key,
        symbols=symbols, natoms=natoms, natypes=natypes)