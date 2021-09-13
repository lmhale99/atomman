# coding: utf-8

def get_crystal_prototypes(self, name=None, id=None, key=None, commonname=None,
                           prototype=None, pearson=None, strukturbericht=None,
                           sg_number=None, sg_hm=None, sg_schoenflies=None,
                           crystalfamily=None, natypes=None, local=None,
                           remote=None, refresh_cache=False, return_df=False,
                           verbose=False):
    """
    Gets all matching crystal prototypes from the database.
    
    Parameters
    ----------
    name : str or list
        The record name(s) to parse by.  For crystal prototype records, the
        names should correspond to the id.
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
    crystalfamily : str, optional
        The crystal system family to limit the search by. 
    natypes : int, optional
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
    """

    return self.get_records(
        style='crystal_prototype', name=name, local=local, remote=remote,
        refresh_cache=refresh_cache, return_df=return_df, verbose=verbose,
        id=id, key=key, commonname=commonname, prototype=prototype,
        pearson=pearson, strukturbericht=strukturbericht,
        sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
        crystalfamily=crystalfamily, natypes=natypes)    

def get_crystal_prototype(self, name=None, id=None, key=None, commonname=None,
                          prototype=None, pearson=None, strukturbericht=None,
                          sg_number=None, sg_hm=None, sg_schoenflies=None,
                          crystalfamily=None, natypes=None, local=None,
                          remote=None, prompt=True, refresh_cache=False,
                          verbose=False):
    
    """
    Gets exactly one matching crystal prototype from the database.
    
    Parameters
    ----------
    name : str or list
        The record name(s) to parse by.  For crystal prototype records, the
        names should correspond to the id.
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
    crystalfamily : str or list, optional
        The crystal system family to limit the search by. 
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
    """
    return self.get_record(
        style='crystal_prototype', name=name, local=True, remote=False,
        prompt=prompt, refresh_cache=refresh_cache, verbose=verbose,
        id=id, key=key, commonname=commonname, prototype=prototype,
        pearson=pearson, strukturbericht=strukturbericht,
        sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
        crystalfamily=crystalfamily, natypes=natypes)

def retrieve_crystal_prototype(self, name=None, dest=None, id=None, key=None,
                               commonname=None, prototype=None, pearson=None,
                               strukturbericht=None, sg_number=None,
                               sg_hm=None, sg_schoenflies=None,
                               crystalfamily=None, natypes=None,
                               local=None, remote=None, prompt=True,
                               format='json', indent=4, refresh_cache=False,
                               verbose=False):
    """
    Gets a single matching crystal prototype from the database and saves it to a
    file based on the record's name.

    Parameters
    ----------
    name : str or list, optional
        The name(s) of records to limit the search by.
    dest : path, optional
        The parent directory where the record will be saved to.  If not given,
        will use the current working directory.
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
    crystalfamily : str or list, optional
        The crystal system family to limit the search by. 
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
        style='crystal_prototype', name=name, dest=dest, local=local,
        remote=remote, prompt=prompt, format=format,
        indent=indent, refresh_cache=refresh_cache, verbose=verbose,
        id=id, key=key, commonname=commonname, prototype=prototype,
        pearson=pearson, strukturbericht=strukturbericht,
        sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
        crystalfamily=crystalfamily, natypes=natypes)

def download_crystal_prototypes(self, name=None, id=None, key=None,
                                commonname=None, prototype=None,
                                pearson=None, strukturbericht=None,
                                sg_number=None, sg_hm=None, sg_schoenflies=None,
                                crystalfamily=None, natypes=None,
                                overwrite=False, return_records=False,
                                verbose=False):
    """
    Downloads crystal prototypes from the remote to the local.
    
    Parameters
    ----------
    name : str or list
        The record name(s) to parse by.  For crystal prototype records, the
        names should correspond to the id.
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
    crystalfamily : str or list, optional
        The crystal system families to limit the search by. 
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
    """

    return self.download_records(
        style='crystal_prototype', name=name, overwrite=overwrite,
        return_records=return_records, verbose=verbose,
        id=id, key=key, commonname=commonname, prototype=prototype,
        pearson=pearson, strukturbericht=strukturbericht,
        sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
        crystalfamily=crystalfamily, natypes=natypes)