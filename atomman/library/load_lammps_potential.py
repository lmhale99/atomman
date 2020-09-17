# coding: utf-8

from . import Database
def load_lammps_potential(id=None, key=None, potid=None, potkey=None,
                   status='active', pair_style=None, elements=None,
                   symbols=None, verbose=False, getfiles=False,
                   database=None,
                   localpath=None, local=True, remote=True):
    """
    Gets a LAMMPS potential from the iprPy library or by downloading from
    potentials.nist.gov if a local copy is not found.  Will issue a prompt if
    multiple LAMMPS potentials match the parameters given.
    
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
    elements : str or list, optional
        The included elemental model(s) to limit the search by.
    symbols : str or list, optional
        The included symbol model(s) to limit the search by.
    verbose: bool, optional
        If True, informative print statements will be used.
    getfiles : bool, optional
        If True, then the parameter files for the matching potentials
        will also be retrieved and copied to the working directory.
        If False (default) and the parameter files are in the library,
        then the returned objects' pot_dir path will be set appropriately.
    database : atomman.library.Database, optional
        A pre-existing Database class to use for accessing the records.  Useful
        if multiple potentials are being retrieved as records can be pre-loaded
        for efficiency.  If not given, a new database object will be
        initialized based on the parameters below.
    localpath : str, optional
        The path to the local database directory to use.  If not given,  will
        use the default potentials database path.  Ignored if database is
        given.
    local : bool, optional
        Indicates if the local version of the library should be searched.
        Default value is True.  Ignored if database is given.
    remote : bool, optional
        Indicates if the remote version of the library should be searched.
        Default value is True.  Ignored if database is given.
        
    Returns
    -------
    potentials.PotentialLAMMPS
        The potential object to use.

    Raises
    ------
    ValueError
        If no matching LAMMPS potentials are found.
    """
    # Create Database object and load if needed
    if database is None:
        if local is True:
            database = Database(load='lammps_potentials', status=status,
                                localpath=localpath, local=local,
                                remote=remote, verbose=verbose)
        else:
            database = Database(local=local, remote=remote, verbose=verbose)
    
    matches = database.get_lammps_potentials(id=id, key=key,
                                             potid=potid, potkey=potkey,
                                             pair_style=pair_style, status=status,
                                             elements=elements, symbols=symbols,
                                             verbose=verbose, getfiles=False)

    if len(matches) == 1:
        potential = matches[0]
    elif len(matches) > 1:
        print('Multiple matching LAMMPS potentials found')
        for i, match in enumerate(matches):
            print(i+1, match.id)
        index = int(input('Please select one:')) - 1
        potential = matches[index]
    else:
        raise ValueError('No matching LAMMPS potentials found')

    if getfiles is True:
        return database.get_lammps_potential(id=potential.id, getfiles=True)
    else:
        return potential