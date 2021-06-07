# coding: utf-8

def get_crystal_prototypes(self, local=None, remote=None, name=None, id=None,
                           key=None, commonname=None, prototype=None,
                           pearson=None, strukturbericht=None, sg_number=None,
                           sg_hm=None, sg_schoenflies=None, crystalfamily=None,
                           natypes=None, return_df=False, verbose=False):
    """
    Get all matching crystal prototypes from the database.
    
    Parameters
    ----------
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
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
    crystalfamily : str, optional
        The crystal system family to limit the search by. 
    natypes : int, optional
        The number(s) of unique atom types to limit the search by.
    return_df : bool, optional
        If True, then the corresponding pandas.Dataframe of metadata
        will also be returned.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    return self.get_records('crystal_prototype', local=local, remote=remote, name=name,
                            id=id, key=key, commonname=commonname, prototype=prototype,
                            pearson=pearson, strukturbericht=strukturbericht,
                            sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
                            crystalfamily=crystalfamily, natypes=natypes,
                            return_df=return_df, verbose=verbose)    

def get_crystal_prototype(self, local=None, remote=None, name=None, id=None,
                          key=None, commonname=None, prototype=None,
                          pearson=None, strukturbericht=None, sg_number=None,
                          sg_hm=None, sg_schoenflies=None, crystalfamily=None,
                          natypes=None, prompt=True, verbose=False):
    
    """
    Retrieves exactly one matching crystal prototype from the database.
    
    Parameters
    ----------
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
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
    crystalfamily : str, optional
        The crystal system family to limit the search by. 
    natypes : int, optional
        The number(s) of unique atom types to limit the search by.
    prompt : bool
        If prompt=True (default) then a screen input will ask for a selection
        if multiple matching potentials are found.  If prompt=False, then an
        error will be thrown if multiple matches are found.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    return self.get_record('crystal_prototype', local=True, remote=False,
                           name=name, id=id, key=key, commonname=commonname,
                           prototype=prototype, pearson=pearson,
                           strukturbericht=strukturbericht, sg_number=sg_number,
                           sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
                           crystalfamily=crystalfamily, natypes=natypes,
                           prompt=prompt, verbose=verbose)
        

def download_crystal_prototypes(self, name=None, id=None, key=None,
                                commonname=None, prototype=None,
                                pearson=None, strukturbericht=None,
                                sg_number=None, sg_hm=None, sg_schoenflies=None,
                                crystalfamily=None, natypes=None,
                                overwrite=False, verbose=False):
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
    crystalfamily : str, optional
        The crystal system family to limit the search by. 
    natypes : int, optional
        The number(s) of unique atom types to limit the search by.
    overwrite : bool, optional
        Flag indicating if any existing local records with names matching
        remote records are updated (True) or left unchanged (False).  Default
        value is False.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    self.download_records('crystal_prototype', name=name,
                          id=id, key=key, commonname=commonname, prototype=prototype,
                          pearson=pearson, strukturbericht=strukturbericht,
                          sg_number=sg_number, sg_hm=sg_hm, sg_schoenflies=sg_schoenflies,
                          crystalfamily=crystalfamily, natypes=natypes,
                          overwrite=overwrite, verbose=verbose)