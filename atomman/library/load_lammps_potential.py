# coding: utf-8

from . import Database
def load_lammps_potential(id=None, key=None, potid=None, potkey=None,
                   status='active', pair_style=None, element=None,
                   symbol=None, verbose=False, get_files=False,
                   database=None,
                   localpath=None, local=True, remote=True):
    """
    Gets a LAMMPS potential from the iprPy library or by downloading from
    potentials.nist.gov if a local copy is not found.  Will raise an error
    if none or multiple matching potentials are found.
    
    Parameters
    ----------
    id : str or list, optional
        The id value(s) to limit the search by.
    key : str or list, optional
        The key value(s) to limit the search by.
    potid : str or list, optional
        The potid value(s) to limit the search by.
    potkey : str or list, optional
        The potkey value(s) to limit the search by.
    status : str or list, optional
        The status value(s) to limit the search by.
    pair_style : str or list, optional
        The pair_style value(s) to limit the search by.
    element : str or list, optional
        The included elemental model(s) to limit the search by.
    symbol : str or list, optional
        The included symbol model(s) to limit the search by.
    verbose: bool, optional
        If True, informative print statements will be used.
    get_files : bool, optional
        If True, then the parameter files for the matching potentials
        will also be retrieved and copied to the working directory.
        If False (default) and the parameter files are in the library,
        then the returned objects' pot_dir path will be set appropriately.
    local : bool, optional
        Indicates if the local version of the library should be searched.
        Default value is True.
    remote : bool, optional
        Indicates if the remote version of the library should be searched.
        Default value is True.
        
    Returns
    -------
    Potential
        The potential object to use.
    """
    # Create Database object and load if needed
    if database is None:
        if local is True:
            database = Database(load='lammps_potentials', localpath=localpath,
                                local=local, remote=remote, verbose=verbose)
        else:
            database = Database(local=local, remote=remote, verbose=verbose)
    
    return database.get_lammps_potential(id=id, key=key, potid=potid, potkey=potkey,
                                 status=status, pair_style=pair_style, element=element,
                                 symbol=symbol, verbose=verbose, get_files=get_files)