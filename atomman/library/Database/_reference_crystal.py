# coding: utf-8

# Standard Python imports
from pathlib import Path
from typing import Optional, Union, Tuple

from yabadaba.record import Record

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# https://docs.python-requests.org/
import requests

# atomman imports
from ... import Atoms, Box, System, load
from ..record import load_record
from ...tools import aslist

def get_reference_crystals(self,
                           name: Union[str, list, None] = None,
                           key: Union[str, list, None] = None,
                           id: Union[str, list, None] = None,
                           sourcename: Union[str, list, None] = None,
                           sourcelink: Union[str, list, None] = None,
                           crystalfamily: Union[str, list, None] = None,
                           composition: Union[str, list, None] = None,
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
    Gets all matching reference crystals from the database.
    
    Parameters
    ----------
    name : str or list, optional
        The record name(s) to parse by.  For reference crystal records, the
        names should correspond to the id.
    id : str or list, optional
        The record id(s) to parse by.  For reference crystal records, the
        id are letters identifying the source database "mp-", "mvc-", or "oqmd-",
        followed by the source database's identification number.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    sourcename : str or list, optional
        The full name of source DFT databases to limit the search by.
        "Materials Project" or "Open Quantum Materials Database".            
    sourcelink : str or list, optional
        The web link of the source DFT databases to limit the search by.
    crystalfamily : str or list, optional
        The crystal system families to limit the search by.
    composition : str or list, optional
        The reduced compositions of the structures to limit the search by.
        Element symbols are sorted alphabetically.
    symbols : str or list, optional
        Element symbols in the crystal to limit the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
    natypes : int or list, optional
        The number(s) of unique atom types to limit the search by.
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
        style='reference_crystal', name=name, local=local, remote=remote, 
        refresh_cache=refresh_cache, return_df=return_df, verbose=verbose,
        key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
        crystalfamily=crystalfamily, composition=composition,
        symbols=symbols, natoms=natoms, natypes=natypes)

def promptfxn(df):
    """
    promptfxn for reference crystals
    """
    header = '#  id           comp     c-family     natoms alat'
    print(header)
    
    js = df.sort_values(['composition', 'natoms', 'crystalfamily', 'a']).index
    for i, j in enumerate(js):
        crystal = df.loc[j]
        row =  f'{i+1:2} {crystal.id:12.12} '
        row += f'{crystal.composition:8.8} '
        row += f'{crystal.crystalfamily:12.12} '
        row += f'{crystal.natoms:5d} '
        row += f'{crystal.a:7.4f}'
        print(row)
    
    i = int(input('Please select one:')) - 1
    if i < 0 or i >= len(js):
        raise ValueError('Invalid selection')

    return js[i]

def get_reference_crystal(self, 
                          name: Union[str, list, None] = None,
                          key: Union[str, list, None] = None,
                          id: Union[str, list, None] = None,
                          sourcename: Union[str, list, None] = None,
                          sourcelink: Union[str, list, None] = None,
                          crystalfamily: Union[str, list, None] = None,
                          composition: Union[str, list, None] = None,
                          symbols: Union[str, list, None] = None,
                          natoms: Union[int, list, None] = None,
                          natypes: Union[int, list, None] = None,
                          local: Optional[bool] = None,
                          remote: Optional[bool] = None,
                          prompt: bool = True,
                          refresh_cache: bool = False,
                          verbose: bool = False) -> Record:
    """
    Gets exactly one matching reference crystal from the database.
    
    Parameters
    ----------
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    name : str or list, optional
        The record name(s) to parse by.  For reference crystal records, the
        names should correspond to the id.
    id : str or list, optional
        The record id(s) to parse by.  For reference crystal records, the
        id are letters identifying the source database "mp-", "mvc-", or "oqmd-",
        followed by the source database's identification number.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    sourcename : str or list, optional
        The full name of source DFT databases to limit the search by.
        "Materials Project" or "Open Quantum Materials Database".            
    sourcelink : str or list, optional
        The web link of the source DFT databases to limit the search by.
    crystalfamily : str or list, optional
        The crystal system families to limit the search by.
    composition : str or list, optional
        The reduced compositions of the structures to limit the search by.
        Element symbols are sorted alphabetically.
    symbols : str or list, optional
        Element symbols in the crystal to limit the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
    natypes : int or list, optional
        The number(s) of unique atom types to limit the search by.
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
        style='reference_crystal', name=name, local=local, remote=remote, prompt=prompt,
        promptfxn=promptfxn, refresh_cache=refresh_cache, verbose=verbose,
        key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
        crystalfamily=crystalfamily, composition=composition,
        symbols=symbols, natoms=natoms, natypes=natypes)

def retrieve_reference_crystal(self,
                               name: Union[str, list, None] = None,
                               dest: Optional[Path] = None,
                               key: Union[str, list, None] = None,
                               id: Union[str, list, None] = None,
                               sourcename: Union[str, list, None] = None,
                               sourcelink: Union[str, list, None] = None,
                               crystalfamily: Union[str, list, None] = None,
                               composition: Union[str, list, None] = None,
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
    id : str or list, optional
        The record id(s) to parse by.  For reference crystal records, the
        id are letters identifying the source database "mp-", "mvc-", or "oqmd-",
        followed by the source database's identification number.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    sourcename : str or list, optional
        The full name of source DFT databases to limit the search by.
        "Materials Project" or "Open Quantum Materials Database".            
    sourcelink : str or list, optional
        The web link of the source DFT databases to limit the search by.
    crystalfamily : str or list, optional
        The crystal system families to limit the search by.
    composition : str or list, optional
        The reduced compositions of the structures to limit the search by.
        Element symbols are sorted alphabetically.
    symbols : str or list, optional
        Element symbols in the crystal to limit the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
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
        style='reference_crystal', name=name, dest=dest, local=local,
        remote=remote, prompt=prompt, promptfxn=promptfxn, format=format,
        indent=indent, refresh_cache=refresh_cache, verbose=verbose,
        key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
        crystalfamily=crystalfamily, composition=composition,
        symbols=symbols, natoms=natoms, natypes=natypes)

def download_reference_crystals(self, 
                                name: Union[str, list, None] = None,
                                key: Union[str, list, None] = None,
                                id: Union[str, list, None] = None,
                                sourcename: Union[str, list, None] = None,
                                sourcelink: Union[str, list, None] = None,
                                crystalfamily: Union[str, list, None] = None,
                                composition: Union[str, list, None] = None,
                                symbols: Union[str, list, None] = None,
                                natoms: Union[int, list, None] = None,
                                natypes: Union[int, list, None] = None,
                                overwrite: bool = False,
                                return_records: bool = False,
                                verbose: bool = False) -> Optional[np.ndarray]:
    """
    Download reference records from the remote and save to localpath.
    
    Parameters
    ----------
    name : str or list, optional
        The record name(s) to parse by.  For reference crystal records, the
        names should correspond to the id.
    id : str or list, optional
        The record id(s) to parse by.  For reference crystal records, the
        id are letters identifying the source database "mp-", "mvc-", or "oqmd-",
        followed by the source database's identification number.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    sourcename : str or list, optional
        The full name of source DFT databases to limit the search by.
        "Materials Project" or "Open Quantum Materials Database".            
    sourcelink : str or list, optional
        The web link of the source DFT databases to limit the search by.
    crystalfamily : str or list, optional
        The crystal system families to limit the search by.
    composition : str or list, optional
        The reduced compositions of the structures to limit the search by.
        Element symbols are sorted alphabetically.
    symbols : str or list, optional
        Element symbols in the crystal to limit the search by.
    natoms : int or list, optional
        The number of unique atoms in the crystal's unit cell to limit the
        search by.
    natypes : int or list, optional
        The number(s) of unique atom types to limit the search by.
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
        key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
        crystalfamily=crystalfamily, composition=composition,
        symbols=symbols, natoms=natoms, natypes=natypes)

def fetch_reference_crystal(self,
                            id: str,
                            api_key: Optional[str] = None,
                            api_key_file: Union[str, Path, None] = None,
                            local: Optional[bool] = None,
                            remote: Optional[bool] = None, 
                            refresh_cache: bool = False,
                            verbose: bool = False) -> Record:
    """
    Retrieves a single reference crystal.  First, the database is checked
    for matches with the DOI, then with the record name.  If no matches are found
    in the database, then the corresponding crystal structure is downloaded from
    the source database.

    Parameters
    ----------
    id : str
        The reference crystal's unique id.  Combines a database tag "mp-" or
        "oqmd-" and the DFT database's entry id.
    api_key : str, optional
        The user's Materials Project API key given as a str.  Either api_key
        or api_key_file are required to retrieve from Materials Project.
    api_key_file : str or Path, optional
        The path to a file containing only the user's Materials Project API
        key.  Either api_key or api_key_file are required to retrieve from
        Materials Project.
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    yabadaba.record.Record
        The matching reference crystal record.
    """
    if local is not False or remote is not False:
        # Try fetching based on doi
        try:
            return self.get_reference_crystal(id=id, local=local, remote=remote,
                                              refresh_cache=refresh_cache, verbose=verbose)
        except:
            pass

    # Fetch from source database
    if 'oqmd-' in id:
        record = self.fetch_oqmd_crystal(id)
        if verbose:
            print('Crystal retrieved from OQMD')
    else:
        record = self.fetch_mp_crystal(id, api_key=api_key, api_key_file=api_key_file)
        if verbose:
            print('Crystal retrieved from Materials Project')

    return record

def fetch_mp_crystals(self, 
                      id: Union[str, list],
                      api_key: Optional[str] = None,
                      api_key_file: Union[Path, str, None] = None) -> list:
    """
    Retrieves reference crystals from Materials Project based on id(s).

    Parameters
    ----------
    id : str or list
        The structure id(s) of the crystals to retrieve.
    api_key : str, optional
        The user's Materials Project API key given as a str.  Either api_key
        or api_key_file are required.
    api_key_file : str or Path, optional
        The path to a file containing only the user's Materials Project API
        key.  Either api_key or api_key_file are required.
    
    Returns
    -------
    list
        All matching reference crystals retrieved.
    """
    def load_pymatgen_structure_json(structure: dict):
    
        # Build box
        a = structure['lattice']['a']
        b = structure['lattice']['b']
        c = structure['lattice']['c']
        alpha = structure['lattice']['alpha']
        beta = structure['lattice']['beta']
        gamma = structure['lattice']['gamma']
        box = Box(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)

        # Init symbols, atype and pos
        symbols = []
        nsites = len(structure['sites'])
        atype = np.empty(nsites, dtype=int)
        pos = np.empty((nsites, 3))

        # Extract symbols, atype and pos from sites
        for i, site in enumerate(structure['sites']):
            if len(site['species']) > 1:
                raise ValueError('Can only do structures where each site is fully occupied by a single element')
            element = site['species'][0]['element']
            if element not in symbols:
                symbols.append(element)
            atype[i] = symbols.index(element) + 1
            pos[i] = np.array(site['abc'])

        # Build atoms and ucell system
        atoms = Atoms(atype=atype, pos=pos)
        ucell = System(atoms=atoms, box=box, symbols=symbols, scale=True)

        return ucell
    
    # Check api_key(_file) values
    if api_key_file is not None:
        if api_key is not None:
            raise ValueError('Cannot give both api_key and api_key_value')
        with open(Path(api_key_file)) as f:
            api_key = f.read().strip()
    elif api_key is None:
        raise ValueError('either api_key or api_key_file must be given')
            
    
    # Build query headers and params
    headers = {"x-api-key": api_key}
    params = {
        "material_ids": ','.join(aslist(id)),
        "deprecated": False,
        '_per_page': 100,
        '_skip': 0,
        '_limit': 100,
        '_fields': 'material_id,symmetry,structure',
        '_all_fields': False,
        'license': 'BY-C'
    }
    
    # Perform REST call
    r = requests.get('https://api.materialsproject.org/materials/core/', headers=headers, params=params)
    data = r.json()['data']
    
    records = []
    for crystal in data:
        mp_id = crystal['material_id']
        
        # Load structure
        ucell = load_pymatgen_structure_json(crystal['structure'])
        
        # Standardize the cell
        ucell = ucell.dump('standardize_cell')
        
        # Build record content
        # Build basic record content
        record = load_record('reference_crystal', name=mp_id)
        record.sourcename = "Materials Project"
        record.sourcelink = "https://materialsproject.org/"
        record.ucell = ucell
        records.append(record)
    
    return records

def fetch_mp_crystal(self,
                     id: str,
                     api_key: Optional[str] = None,
                     api_key_file: Union[Path, str, None] = None) -> Record:
    """
    Retrieves a single reference crystal from Materials Project based on id.

    Parameters
    ----------
    id : str
        The structure id of the crystal to retrieve.
    api_key : str, optional
        The user's Materials Project API key given as a str.  Either api_key
        or api_key_file are required.
    api_key_file : str or Path, optional
        The path to a file containing only the user's Materials Project API
        key.  Either api_key or api_key_file are required.
    
    Returns
    -------
    Record
        The retrieved reference crystal.
    """
    # Read api_key from a file
    if Path(api_key).is_file():
        with open(api_key) as f:
            api_key = f.read().strip()

    records = self.fetch_mp_crystals(id, api_key=api_key, api_key_file=api_key_file)
    if len(records) == 1:
        return records[0]
    else:
        raise ValueError('Exactly one record not found?')

def fetch_oqmd_crystal(self, id: str) -> Record:
    """
    Retrieves a single reference crystal from OQMD based on id.

    Parameters
    ----------
    id : str
        The OQMD entry number with "oqmd-" prefix.

    Returns
    -------
    Record
        The retrieved reference crystal.
    """
    
    # Build basic record content
    record = load_record('reference_crystal', name=id)
    record.sourcename = "Open Quantum Materials Database"
    record.sourcelink = "http://oqmd.org/"
    
    # Parse entry number page for structure number
    entry_number = id.replace('oqmd-', '')
    entry_r = requests.get(f'http://oqmd.org/materials/entry/{entry_number}')
    entry_html = entry_r.text        
    start = entry_html.index('href="/materials/structure/') + len('href="/materials/structure/')
    end = start + entry_html[start:].index('"')
    structure_number = entry_html[start:end]

    # Try retrieving poscar of conventional then primitive cells
    try:
        structure_url = f'http://oqmd.org/materials/export/conventional/poscar/{structure_number}'
        structure_r = requests.get(structure_url)
        structure_r.raise_for_status()
    except:
        try:
            structure_url = f'http://oqmd.org/materials/export/primitive/poscar/{structure_number}' 
            structure_r = requests.get(structure_url)
            structure_r.raise_for_status()
        except:
            raise ValueError('Failed to find the poscar file for the structure')
            
    # Load ucell
    record.ucell = load('poscar', structure_r.text).normalize()
    
    return record

def save_reference_crystal(self,
                           crystal: Record,
                           overwrite: bool = False,
                           verbose: bool = False):
    """
    Saves a reference crystal to the local database.
    
    Parameters
    ----------
    crystal : ReferenceCrystal
        The record to save.  
    overwrite : bool, optional
        Indicates what to do when a matching record is found in the local
        location.  If False (default), then the record is not updated.  If
        True, then the record is updated.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    self.save_record(record=crystal, overwrite=overwrite, verbose=verbose)