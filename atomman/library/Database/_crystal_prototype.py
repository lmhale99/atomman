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

from .. import CrystalPrototype
from ...tools import aslist, screen_input

@property
def crystal_prototypes(self):
    """list or None: Loaded CrystalPrototype objects"""
    return self.__crystal_prototypes

@property
def crystal_prototypes_df(self):
    """pandas.DataFrame or None: Metadata for loaded CrystalPrototype objects"""
    return self.__crystal_prototypes_df

def _no_load_crystal_prototypes(self):
    """Initializes properties if load_crystal_prototypes is not called"""
    self.__crystal_prototypes = None
    self.__crystal_prototypes_df = None

def load_crystal_prototypes(self, localpath=None, local=None, remote=None, verbose=False):
    """
    Loads crystal prototypes from the database, first checking localpath, then
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
    prototypes = {}
    
    # Set localpath, local, and remote as given here or during init
    if localpath is None:
        localpath = self.localpath
    if local is None:
        local = self.local
    if remote is None:
        remote = self.remote
    
    # Check localpath first
    if local is True and localpath is not None:
        for fname in Path(localpath, 'crystal_prototype').glob('*'):
            if fname.suffix in ['.xml', '.json']:
                
                with open(fname, encoding='UTF-8') as f:
                    proto = CrystalPrototype(model=f.read())
                prototypes[proto.id] = proto
                
        if verbose:
            print(f'Loaded {len(prototypes)} local crystal prototypes')
    numlocal = len(prototypes)

    # Load remote
    if remote is True:
        try:
            records = self.cdcs.query(template='crystal_prototype')
        except:
            if verbose:
                print('Failed to load crystal prototypes from remote')
        else:
            if verbose:
                print(f'Loaded {len(records)} remote crystal prototypes')
            
            for i in range(len(records)):
                record = records.iloc[i]
                proto = CrystalPrototype(record.xml_content)
                if proto.id not in prototypes:
                    prototypes[proto.id] = proto

            if verbose and len(prototypes) > 0:
                print(f' - {len(prototypes) - numlocal} new')
    
    # Build crystal_prototypes and crystal_prototypes_df
    if len(prototypes) > 0:
        protos = np.array(list(prototypes.values()))
        protodicts = []
        for proto in prototypes.values():
            protodicts.append(proto.asdict())
        
        self.__crystal_prototypes_df = pd.DataFrame(protodicts).sort_values('id')
        self.__crystal_prototypes = protos[self.crystal_prototypes_df.index]
        self.__crystal_prototypes_df.reset_index(drop=True)
    else:
        if verbose:
            print('No crystal prototypes loaded')
        self.__crystal_prototypes = None
        self.__crystal_prototypes_df = None
        
def get_crystal_prototypes(self, id=None, key=None, name=None, prototype=None,
                   pearson=None, strukturbericht=None, sg_number=None,
                   sg_hm=None, sg_schoenflies=None, keyword=None,
                   verbose=False):
    """
    Get all matching crystal prototypes from the database.
    
    Parameters
    ----------
    id : str or list, optional
        Prototype ID(s) to search for.  These are unique identifiers for each
        prototype based on comm.
    key : str or list, optional
        UUID4 key(s) to search for.  Each entry has a unique random-generated
        UUID4 key.
    name : str or list, optional
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
    keyword : str, optional
        If given, will limit the search to all records that contain the keyword
        substring.  Cannot be combined with any of the above parameters.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    list of Potential
        The matching potentials
    """

    if keyword is not None:
        try:
            assert id is None
            assert key is None
            assert name is None
            assert prototype is None
            assert pearson is None
            assert strukturbericht is None
            assert sg_hm is None
            assert sg_number is None
            assert sg_schoenflies is None
        except:
            raise ValueError('keyword cannot be combined with the other search limiting paramters')

    # Check loaded values if available
    if self.crystal_prototypes_df is not None:
        
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
            
            prototypes = self.crystal_prototypes[self.crystal_prototypes_df.apply(strmatch, args=['id', id], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['key', key], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['name', name], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['prototype', prototype], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['pearson', pearson], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['strukturbericht', strukturbericht], axis=1)
                                                &self.crystal_prototypes_df.apply(intmatch, args=['sg_number', sg_number], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['sg_hm', sg_hm], axis=1)
                                                &self.crystal_prototypes_df.apply(strmatch, args=['sg_schoenflies', sg_schoenflies], axis=1)]
        else:
            
            prototypes = []
            for proto in self.crystal_prototypes:
                if keyword.lower() in proto.asmodel().json().lower():
                    prototypes.append(proto)

        if verbose:
            print(len(prototypes), 'matching crystal prototypes found from loaded records')
        return prototypes

    # Check remote values if no loaded values
    else:
        if keyword is None:
            # Build Mongoquery
            mquery = {}

            # Add id query
            if id is not None:
                id = aslist(id)
                mquery['crystal-prototype.id'] = {'$in': id}

            # Add key query
            if key is not None:
                key = aslist(key)
                mquery['crystal-prototype.key'] = {'$in': key}

            # Add name query
            if name is not None:
                name = aslist(name)
                mquery['crystal-prototype.name'] = {'$in': name}

            # Add prototype query
            if prototype is not None:
                prototype = aslist(prototype)
                mquery['crystal-prototype.prototype'] = {'$in': prototype}

            # Add pearson query
            if pearson is not None:
                pearson = aslist(pearson)
                mquery['crystal-prototype.Pearson-symbol'] = {'$in': pearson}

            # Add strukturbericht query
            if strukturbericht is not None:
                strukturbericht = aslist(strukturbericht)
                mquery['crystal-prototype.Strukturbericht'] = {'$in': strukturbericht}

            # Add sg_number query
            if sg_number is not None:
                sg_number = aslist(sg_number)
                for i in range(len(sg_number)):
                    sg_number[i] = int(sg_number[i])
                mquery['crystal-prototype.space-group.number'] = {'$in': sg_number}
            
            # Add sg_hm query
            if sg_hm is not None:
                sg_hm = aslist(sg_hm)
                mquery['crystal-prototype.space-group.Hermann-Maguin'] = {'$in': sg_hm}

            # Add sg_schoenflies query
            if sg_schoenflies is not None:
                sg_schoenflies = aslist(sg_schoenflies)
                mquery['crystal-prototype.space-group.Schoenflies'] = {'$in': sg_schoenflies}

            matches = self.cdcs.query(template='crystal_prototype', mongoquery=mquery)

        else:
            matches = self.cdcs.query(template='crystal_prototype', keyword=keyword)

        if verbose:
            print(len(matches), 'matching crystal prototypes found from remote database')

        if len(matches) > 0:
            matches = matches.sort_values('id').reset_index(drop=True)

            def makeprototypes(series):
                return CrystalPrototype(model=series.xml_content)

            return matches.apply(makeprototypes, axis=1).values
        else:
            return np.array([])

def get_crystal_prototype(self, id=None, key=None, name=None, prototype=None,
                   pearson=None, strukturbericht=None, sg_number=None,
                   sg_hm=None, sg_schoenflies=None, keyword=None,
                   verbose=False):
    
    prototypes = self.get_crystal_prototypes(id=id, key=key, name=name, prototype=prototype,
                   pearson=pearson, strukturbericht=strukturbericht, sg_number=sg_number,
                   sg_hm=sg_hm, sg_schoenflies=sg_schoenflies, keyword=keyword,
                   verbose=verbose)

    # Check number of matches and select
    if len(prototypes) == 1:
        return prototypes[0]
    
    elif len(prototypes) > 1:
        print('Multiple matching crystal prototypes found.')
        for i, prototype in enumerate(prototypes):
            print(i+1, prototype.id)
        choice = int(screen_input('Select which one:'))
        return prototypes[choice-1]
    
    else:
        raise ValueError('No matching prototypes found')

def download_crystal_prototypes(self, localpath=None, prototypes=None, format='json',
                        indent=None, verbose=False):
    """
    Download citation records from the remote and save to localpath.
    
    Parameters
    ----------
    localpath : path-like object, optional
        Path to a local directory where the files will be saved to.  If not
        given, will use the localpath value set during object initialization.
    prototypes : list of CrystalPrototype, optional
        A list of crystal prototypes to download. If not given, all prototypes will
        be downloaded.
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

    template = 'crystal_prototype'

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
    if prototypes is None:
        self.download_records(template=template, localpath=localpath,
                              format=format, indent=indent, verbose=verbose)
        
    else:
        prototypes = aslist(prototypes)

        for prototype in prototypes:
            fname = Path(save_directory, f'{prototype.id}.{format}')
            if format == 'xml':
                with open(fname, 'w', encoding='UTF-8') as f:
                    prototype.asmodel().xml(fp=f, indent=indent)
            elif format == 'json':
                with open(fname, 'w', encoding='UTF-8') as f:
                    prototype.asmodel().json(fp=f, indent=indent)
        if verbose:
            print(f'Copied {len(prototypes)} records to localpath')

def upload_crystal_prototype(self, prototype, workspace=None, verbose=False):
    raise NotImplementedError('To be done...')