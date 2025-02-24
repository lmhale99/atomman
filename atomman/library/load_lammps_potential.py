# coding: utf-8

# Standard Python imports
from typing import Optional, Union

# atomman imports
from . import Database

def load_lammps_potential(name: Union[str, list, None] = None,
                          key: Union[str, list, None] = None,
                          id: Union[str, list, None] = None,
                          potid: Union[str, list, None] = None,
                          potkey: Union[str, list, None] = None,
                          units: Union[str, list, None] = None,
                          atom_style: Union[str, list, None] = None,
                          pair_style: Union[str, list, None] = None,
                          status: Union[str, list, None] = None,
                          symbols: Union[str, list, None] = None,
                          elements: Union[str, list, None] = None,
                          pot_dir_style: Optional[str] = None,
                          kim_models: Optional[list] = None,
                          kim_api_directory: Optional[str] = None,
                          kim_models_file: Optional[str] = None, 
                          local: Optional[bool] = None,
                          remote: Optional[bool] = None,
                          database: Optional[Database] = None,
                          getfiles: bool = False,
                          prompt: bool = True,
                          verbose: bool = False):
    """
    Loads a LAMMPS potential from the NIST Interatomic Potentials Repository
    or from a local copy of the repository.  Will issue a prompt if
    multiple LAMMPS potentials match the parameters given.
    
    Parameters
    ----------
    name : str or list
        The record name(s) to parse by.  For potential records, the names
        should correspond to the id with a prefix of "potentials." added to it.
    key : str or list
        The unique UUID4 record key(s) to parse by.
    id : str or list
        The unique record id(s) labeling the records to parse by.
    potid : str or list
        The unique UUID4 record key(s) for the associated potential records to
        parse by.
    potkey : str or list
        The unique record id(s) labeling the associated potential records to
        parse by.
    units : str or list
        LAMMPS units option(s) to parse by.
    atom_style : str or list
        LAMMPS pair_style(s) to parse by.
    pair_style : str or list
        LAMMPS pair_style(s) to parse by.
    status : None, str or list
        Limits the search by the status of the LAMMPS implementations:
        "active", "superseded" and/or "retracted".  By default, only active
        implementations are returned.  Giving a value of None will return
        implementations of all statuses.
    symbols : str or list
        Model symbol(s) to parse by.  Typically correspond to elements for
        atomic potential models.
    elements : str or list
        Element(s) in the model to parse by.
    pot_dir_style : str, optional
        Specifies how the pot_dir values will be set for the retrieved LAMMPS
        potentials.  Allowed values are 'working', 'id', and 'local'.
        'working' will set all pot_dir = '', meaning parameter files
        are expected in the working directory when the potential is accessed.
        'id' sets the pot_dir values to match the potential's id.
        'local' sets the pot_dir values to the corresponding local database
        paths where the files are expected to be found.  Default value is
        controlled by settings.
    kim_models : list
        A list of full KIM model ids to build LAMMPS potentials for.
    kim_api_directory : str
        The path to the directory containing a kim-api-collections-management
        executable to use to identify which KIM models are installed.
    kim_models_file : str
        The path to a file containing a list of full KIM model ids to build
        LAMMPS potentials for.
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    database : potentials.Database, optional
        Allows for a previously defined Database object to be used to find
        the potential.  If not given, a new Database object will be used with
        the default local and remote interaction settings.
    getfiles : bool, optional
        If True, then the parameter files for the matching potentials
        will also be copied/downloaded to the potential directory.
    prompt : bool
        If prompt=True (default) then a screen input will ask for a selection
        if multiple matching potentials are found.  If prompt=False, then an
        error will be thrown if multiple matches are found.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    
    Returns
    -------
    potentials.record.PotentialLAMMPS
        The potential object to use.

    Raises
    ------
    ValueError
        If no matching LAMMPS potentials are found.
    """
    # Create Database object and load if needed
    if database is None:
        database = Database()
    
    lmppot = database.get_lammps_potential(name=name, key=key, id=id,
                                           potid=potid, potkey=potkey,
                                           units=units, atom_style=atom_style,
                                           pair_style=pair_style, status=status,
                                           symbols=symbols, elements=elements,
                                           pot_dir_style=pot_dir_style,
                                           kim_models=kim_models,
                                           kim_api_directory=kim_api_directory,
                                           kim_models_file=kim_models_file,
                                           local=local, remote=remote,
                                           prompt=prompt, verbose=verbose)

    if getfiles is True:
        database.get_lammps_potential_files(lmppot, local=local, remote=remote,
                                            verbose=verbose)
    return lmppot
