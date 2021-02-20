# coding: utf-8
# Standard Python libraries
from pathlib import Path


import pandas as pd


import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/atomman
import atomman as am

import potentials

from .. import RelaxedCrystal, CrystalPrototype
from ...tools import aslist, screen_input

@property
def relaxed_crystals(self):
    """list or None: Loaded RelaxedCrystal objects"""
    return self.__relaxed_crystals

@property
def relaxed_crystals_df(self):
    """pandas.DataFrame or None: Metadata for loaded RelaxedCrystal objects"""
    return self.__relaxed_crystals_df

def _no_load_relaxed_crystals(self):
    """Initializes properties if load_relaxed_crystals is not called"""
    self.__relaxed_crystals = None
    self.__relaxed_crystals_df = None

def load_relaxed_crystals(self, localpath=None, local=None, remote=None, verbose=False):
    """
    Loads relaxed crystals from the database, first checking localpath, then
    trying to download from host.
    
    Parameters
    ----------
    localpath : str, optional
        Path to a local directory to check for records first.  If not given,
        will check localpath value set during object initialization.  If not
        given or set during initialization, then only the remote database will
        be loaded.
    local : bool, optional
        Indicates if records in localpath are to be loaded.  If not given,
        will use the local value set during initialization.
    remote : bool, optional
        Indicates if the records in the remote database are to be loaded.
        Setting this to be False is useful/faster if a local copy of the
        database exists.  If not given, will use the local value set during
        initialization.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    crystals = {}
    
    # Set localpath, local, and remote as given here or during init
    if localpath is None:
        localpath = self.localpath
    if local is None:
        local = self.local
    if remote is None:
        remote = self.remote
    
    # Check localpath first
    if local is True and localpath is not None:
        for fname in Path(localpath, 'relaxed_crystal').glob('*'):
            if fname.suffix in ['.xml', '.json']:
                
                with open(fname, encoding='UTF-8') as f:
                    crystal = RelaxedCrystal(model=f.read())
                crystals[crystal.key] = crystal
                
        if verbose:
            print(f'Loaded {len(crystals)} local relaxed crystals')
    numlocal = len(crystals)

    # Load remote
    if remote is True:
        try:
            records = self.cdcs.query(template='relaxed_crystal')
        except:
            if verbose:
                print('Failed to load relaxed crystals from remote')
        else:
            if verbose:
                print(f'Loaded {len(records)} remote relaxed crystals')
            
            for i in range(len(records)):
                record = records.iloc[i]
                crystal = RelaxedCrystal(record.xml_content)
                if crystal.key not in crystals:
                    crystals[crystal.key] = crystal

            if verbose and len(crystals) > 0 and local:
                print(f' - {len(crystals) - numlocal} new')
    
    # Build relaxed_crystals and relaxed_crystals_df
    if len(crystals) > 0:
        crys = np.array(list(crystals.values()))
        crysdicts = []
        for cry in crystals.values():
            crysdicts.append(cry.asdict())
        
        self.__relaxed_crystals_df = pd.DataFrame(crysdicts).sort_values('key')
        self.__relaxed_crystals = crys[self.relaxed_crystals_df.index]
        self.__relaxed_crystals_df.reset_index(drop=True)
    else:
        if verbose:
            print('No relaxed crystals loaded')
        self.__relaxed_crystals = None
        self.__relaxed_crystals_df = None
        
def get_relaxed_crystals(self, key=None, method='dynamic', standing='good',
                   family=None, parent_key=None, potential=None,
                   potential_LAMMPS_id=None, potential_LAMMPS_key=None, potential_id=None, potential_key=None,
                   symbols=None, natoms=None, natypes=None, keyword=None,
                   verbose=False):
    """
    Get all matching relaxed crystals from the database.
    
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    list of RelaxedCrystal
        The matching relaxed crystal records.
    """

    if keyword is not None:
        try:
            assert key is None
            #assert method is None
            #assert standing is None
            assert parent_key is None
            assert family is None
            assert potential is None
            assert potential_LAMMPS_id is None
            assert potential_LAMMPS_key is None
            assert potential_id is None
            assert potential_key is None
            assert symbols is None
            assert natypes is None
            assert natoms is None
        except:
            raise ValueError('keyword cannot be combined with the other search limiting paramters')

    else:
        # Extract potential ids from potential objects
        if potential is not None:
            try:
                assert potential_LAMMPS_id is None
                assert potential_LAMMPS_key is None
                assert potential_id is None
                assert potential_key is None
            except:
                raise ValueError('potential cannot be given with the other potential_* parameters')
            potential_LAMMPS_key = []
            for pot in aslist(potential):
                potential_LAMMPS_key.append(pot.key)

        # Extract family names from prototypes
        if family is not None:
            family = aslist(family)
            for i in range(len(family)):
                if isinstance(family[i], CrystalPrototype):
                    family[i] = family[i].id

    # Check loaded values if available
    if self.relaxed_crystals_df is not None:
        
        if keyword is None:

            def strmatch(series, name, val):
                if val is None:
                    return True
                else:
                    return series[name] in aslist(val)

            def intmatch(series, name, val):
                if val is None:
                    return True
                else:
                    val = aslist(val)
                    for i in range(len(val)):
                        val[i] = int(val[i])
                    return series[name] in val
            
            def listmatch(series, name, val):
                if val is None:
                    return True
                
                elif isinstance(series[name], (list, tuple)):
                    for v in aslist(val):
                        if v not in series[name]:
                            return False
                    return True
                else:
                    return False
            
            crystals = self.relaxed_crystals[self.relaxed_crystals_df.apply(strmatch, args=['key', key], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['method', method], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['standing', standing], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['parent_key', parent_key], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['family', family], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['potential_LAMMPS_id', potential_LAMMPS_id], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['potential_LAMMPS_key', potential_LAMMPS_key], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['potential_id', potential_id], axis=1)
                                             &self.relaxed_crystals_df.apply(strmatch, args=['potential_key', potential_key], axis=1)
                                             &self.relaxed_crystals_df.apply(intmatch, args=['natoms', natoms], axis=1)
                                             &self.relaxed_crystals_df.apply(intmatch, args=['natypes', natypes], axis=1)
                                             &self.relaxed_crystals_df.apply(listmatch, args=['symbols', symbols], axis=1)]
        else:
            
            crystals = []
            for crystal in self.relaxed_crystals:
                if keyword.lower() in crystal.asmodel().json().lower():
                    crystals.append(crystal)

        if verbose:
            print(len(crystals), 'matching relaxed crystals found from loaded records')
        return crystals

    # Check remote values if no loaded values
    else:
        if keyword is None:
            # Build Mongoquery
            mquery = {}

            # Add key query
            if key is not None:
                key = aslist(key)
                mquery['relaxed-crystal.key'] = {'$in': key}

            # Add method query
            if method is not None:
                method = aslist(method)
                mquery['relaxed-crystal.method'] = {'$in': method}

            # Add standing query
            if standing is not None:
                standing = aslist(standing)
                mquery['relaxed-crystal.standing'] = {'$in': standing}

            # Add family query
            if family is not None:
                family = aslist(family)
                mquery['relaxed-crystal.system-info.family'] = {'$in': family}

            # Add parent_key query
            if parent_key is not None:
                parent_key = aslist(parent_key)
                mquery['relaxed-crystal.system-info.parent_key'] = {'$in': parent_key}
            
            # Add potential_LAMMPS_id query
            if potential_LAMMPS_id is not None:
                potential_LAMMPS_id = aslist(potential_LAMMPS_id)
                mquery['relaxed-crystal.potential-LAMMPS.id'] = {'$in': potential_LAMMPS_id}

            # Add potential_LAMMPS_key query
            if potential_LAMMPS_key is not None:
                potential_LAMMPS_key = aslist(potential_LAMMPS_key)
                mquery['relaxed-crystal.potential-LAMMPS.key'] = {'$in': potential_LAMMPS_key}

            # Add potential_id query
            if potential_id is not None:
                potential_id = aslist(potential_id)
                mquery['relaxed-crystal.potential-LAMMPS.potential.id'] = {'$in': potential_id}

            # Add potential_key query
            if potential_key is not None:
                potential_key = aslist(potential_key)
                mquery['relaxed-crystal.potential-LAMMPS.potential.key'] = {'$in': potential_key}

            # Add natoms query
            if natoms is not None:
                natoms = aslist(natoms)
                for i in range(len(natoms)):
                    natoms[i] = int(natoms[i])
                mquery['relaxed-crystal.atomic-system.atoms.natoms'] = {'$in': natoms}

            # Add natypes query
            if natypes is not None:
                natypes = aslist(natypes)
                mquery['$or'] = []
                for natype in natypes:
                    if natype == 1:
                        mquery['$or'].append({"relaxed-crystal.atomic-system.atom-type-symbol":{"$not":{"$type":"array"}}})
                    else:
                        mquery['$or'].append({"relaxed-crystal.atomic-system.atom-type-symbol":{"$size":natype}})

            # Add symbol query
            if symbols is not None:
                symbols = aslist(symbols)
                mquery['relaxed-crystal.atomic-system.atom-type-symbol'] = {'$all': symbols}
            
            matches = self.cdcs.query(template='relaxed_crystal', mongoquery=mquery)

        else:
            matches = self.cdcs.query(template='relaxed_crystal', keyword=keyword)

        if verbose:
            print(len(matches), 'matching relaxed crystals found from remote database')

        if len(matches) > 0:
            def makecrystals(series):
                return RelaxedCrystal(model=series.xml_content)
            crystals = matches.apply(makecrystals, axis=1).values

            crystals_df = []
            for crystal in crystals:
                crystals_df.append(crystal.asdict())
            crystals_df = pd.DataFrame(crystals_df)

            return crystals[crystals_df.sort_values('cohesive_energy').index]
            
        else:
            return np.array([])

def get_relaxed_crystal(self, key=None, method='dynamic', standing='good',
                   family=None, parent_key=None, potential=None,
                   potential_LAMMPS_id=None, potential_LAMMPS_key=None, potential_id=None, potential_key=None,
                   symbols=None, natoms=None, natypes=None, keyword=None,
                   verbose=False):
    """
    Gets a single matching relaxed crystal from the database. If multiple
    matches are found, a selection menu will appear.
    
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    RelaxedCrystal
        A selected relaxed crystal record.

    Raises
    ------
    ValueError
        If no matching records are identified.
    """
    crystals = self.get_relaxed_crystals(key=key, method=method, standing=standing,
                   family=family, parent_key=parent_key, potential=potential,
                   potential_LAMMPS_id=potential_LAMMPS_id, potential_LAMMPS_key=potential_LAMMPS_key,
                   potential_id=potential_id, potential_key=potential_key,
                   symbols=symbols, natoms=natoms, natypes=natypes, keyword=keyword,
                   verbose=verbose)

    # Check number of matches and select
    if len(crystals) == 1:
        return crystals[0]
    
    elif len(crystals) > 1:
        print('Multiple matching relaxed crystals found.')

        # Create list header
        header = '#  family               symbols  alat    Ecoh    '
        if method is None:
            header += 'method  '
        if standing is None:
            header += 'standing'
        print(header)
        
        # Generate options list
        for i in range(len(crystals)):
            crystal = crystals[i]
            
            row = f'{i+1:2} {crystal.family:20.20} {"".join(crystal.ucell.symbols):8.8} {crystal.ucell.box.a:7.4f} {crystal.cohesive_energy:7.4f} '
            
            if method is None:
                row += f'{crystal.method:7.7} '
            if standing is None:
                row += f'{crystal.standing:4.4}'
            print(row)

        choice = int(screen_input('Select which one:'))
        return crystals[choice-1]
    
    else:
        raise ValueError('No matching relaxed crystals found')

def download_relaxed_crystals(self, localpath=None, crystals=None, standing='good',
                              format='json', indent=None, verbose=False):
    """
    Download citation records from the remote and save to localpath.
    
    Parameters
    ----------
    localpath : path-like object, optional
        Path to a local directory where the files will be saved to.  If not
        given, will use the localpath value set during object initialization.
    crystals : list of RelaxedCrystals, optional
        A list of relaxed crystals to download. If not given, all (good)
        crystals will be downloaded.
    standing : str or None, optional
        By default, only records with standing == good will be downloaded.
        Setting this to None will download all crystals.
    format : str, optional
        The file format to save the record files as.  Allowed values are
        'xml' and 'json' (default).
    indent : int, optional
        The indentation spacing size to use for the locally saved record files.
        If not given, the JSON/XML content will be compact.  Ignored if format
        is 'bib'.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    template = 'relaxed_crystal'

    # Handle localpath value
    if localpath is None:
        localpath = self.localpath
    if localpath is None:
        raise ValueError('No local path set to save files to')
    
    # Check format value
    format = format.lower()
    allowed_formats = ['xml', 'json']
    if format not in allowed_formats:
        raise ValueError("Format must be 'xml' or 'json'")

    # Create save directory if needed
    save_directory = Path(localpath, template)
    if not save_directory.is_dir():
        save_directory.mkdir(parents=True)

    for fmt in allowed_formats:
        if fmt != format:
            numexisting = len([fname for fname in save_directory.glob(f'*.{fmt}')])
            if numexisting > 0:
                raise ValueError(f'{numexisting} records of format {fmt} already saved locally')

    # Download if needed
    if crystals is None:
        if standing is not None:
            mquery = {}
            mquery['relaxed-crystal.standing'] = {'$in': standing}
        else:
            mquery = None

        self.download_records(template=template, localpath=localpath, mongoquery=mquery,
                              format=format, indent=indent, verbose=verbose)
        
    else:
        crystals = aslist(crystals)

        for crystal in crystals:
            fname = Path(save_directory, f'{crystal.id}.{format}')
            if format == 'xml':
                with open(fname, 'w', encoding='UTF-8') as f:
                    crystal.asmodel().xml(fp=f, indent=indent)
            elif format == 'json':
                with open(fname, 'w', encoding='UTF-8') as f:
                    crystal.asmodel().json(fp=f, indent=indent)
        if verbose:
            print(f'Copied {len(crystals)} records to localpath')

def upload_relaxed_crystal(self, crystal, workspace=None, verbose=False):
    raise NotImplementedError('To be done...')