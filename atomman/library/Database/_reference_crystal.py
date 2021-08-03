# coding: utf-8

import requests

from ... import load
from ..record import ReferenceCrystal
from ...tools import aslist

def get_reference_crystals(self, local=None, remote=None, name=None, key=None,
                           id=None, sourcename=None, sourcelink=None,
                           crystalfamily=None, composition=None,
                           symbols=None, natoms=None, natypes=None,
                           return_df=False, verbose=False):
    """
    Get all matching reference crystals from the database.
    
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
    return_df : bool, optional
        If True, then the corresponding pandas.Dataframe of metadata
        will also be returned.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    return self.get_records('reference_crystal', local=local, remote=remote, name=name,
                            key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
                            crystalfamily=crystalfamily, composition=composition,
                            symbols=symbols, natoms=natoms, natypes=natypes,
                            return_df=return_df, verbose=verbose)

def get_reference_crystal(self, local=None, remote=None, name=None,
                        key=None, id=None, sourcename=None, sourcelink=None,
                        crystalfamily=None, composition=None,
                        symbols=None, natoms=None, natypes=None,
                        prompt=True, verbose=False):
    """
    Retrieves exactly one matching reference crystal from the database.
    
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    def promptfxn(df):
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

    return self.get_record('reference_crystal', local=local, remote=remote, name=name,
                           key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
                           crystalfamily=crystalfamily, composition=composition,
                           symbols=symbols, natoms=natoms, natypes=natypes,
                           prompt=prompt, promptfxn=promptfxn, verbose=verbose)

def download_reference_crystals(self, name=None, key=None, id=None,
                                sourcename=None, sourcelink=None,
                                crystalfamily=None, composition=None,
                                symbols=None, natoms=None, natypes=None, keyword=None,
                                overwrite=False, verbose=False):
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    self.download_records('relaxed_crystal', name=name,
                          key=key, id=id, sourcename=sourcename, sourcelink=sourcelink,
                          crystalfamily=crystalfamily, composition=composition,
                          symbols=symbols, natoms=natoms, natypes=natypes,
                          overwrite=overwrite, verbose=verbose)

def fetch_reference_crystal(self, id, api_key=None, local=None, remote=None, verbose=False):
    """
    Retrieves a single reference crystal.  First, the database is checked
    for matches with the DOI, then with the record name.  If no matches are found
    in the database, then the corresponding crystal structure is downloaded from
    the source database.

    Parameters
    ----------
    id : str
        The reference crystal's unique id.
    api_key : str, optional
        The user's Materials Project API key. If not given, will use "MAPI_KEY"
        environment variable to fetch records from Materials Project if needed.
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    if local is not False or remote is not False:
        # Try fetching based on doi
        try:
            return self.get_referene_crystal(id=id, local=local, remote=remote, verbose=verbose)
        except:
            pass

    # Fetch from source database
    if 'oqmd-' in id:
        record = self.fetch_oqmd_crystal(id)
        if verbose:
            print(f'Crystal retrieved from OQMD')
    else:
        record = self.fetch_mp_crystal(id, api_key=api_key)
        if verbose:
            print(f'Crystal retrieved from Materials Project')

    return record

def fetch_mp_crystals(self, id, api_key=None):
    """
    Retrieves reference crystals from Materials Project based on id(s).

    Parameters
    ----------
    id : str or list
        The structure id(s) of the crystals to retrieve.
    api_key : str, optional
        The user's Materials Project API key. If not given, will use "MAPI_KEY"
        environment variable.
    """

    # Function-specific imports
    import pymatgen as pmg
    from pymatgen.ext.matproj import MPRester

    # Open connection to Materials Project
    records = []
    with MPRester(api_key) as m:

        # Download missing entries
        try:
            entries = m.query({"material_id": {"$in": aslist(id)}}, ['material_id', 'cif'])
        except:
            raise ValueError('Failed to find Materials Project information')
        else:
            # Convert cif to model and save
            for entry in entries:
                entry_id = entry['material_id']
                struct = pmg.Structure.from_str(entry['cif'], fmt='cif')
                struct = pmg.symmetry.analyzer.SpacegroupAnalyzer(struct).get_conventional_standard_structure()

                # Build record content
                # Build basic record content
                record = ReferenceCrystal(name=entry_id)
                record.sourcename = "Materials Project"
                record.sourcelink = "https://materialsproject.org/"
                record.ucell = load('pymatgen_Structure', struct).normalize()
                records.append(record)
    
    return records

def fetch_mp_crystal(self, id, api_key=None):
    """
    Retrieves a single reference crystal from Materials Project based on id.

    Parameters
    ----------
    id : str
        The structure id of the crystal to retrieve.
    api_key : str, optional
        The user's Materials Project API key. If not given, will use "MAPI_KEY"
        environment variable.
    """
    records = fetch_mp_crystals(id, api_key=api_key)
    if len(records) == 1:
        return records[0]
    else:
        raise ValueError('Exactly one record not found?')

def fetch_oqmd_crystal(self, id):
    """
    Retrieves a single reference crystal from OQMD based on id.

    Parameters
    ----------
    id : str
        The OQMD entry number with "oqmd-" prefix.
    """
    
    # Build basic record content
    record = ReferenceCrystal(name=id)
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

def save_reference_crystal(self, crystal, overwrite=False, verbose=False):
    """
    Saves a citation to the local database.
    
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