# coding: utf-8
        
def get_relaxed_crystals(self, local=None, remote=None, name=None, key=None,
                         method='dynamic', standing='good',
                         family=None, parent_key=None, 
                         potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                         potential_id=None, potential_key=None,
                         symbols=None, natoms=None, natypes=None,
                         return_df=False, verbose=False):
    """
    Get all matching relaxed crystals from the database.
    
    Parameters
    ----------
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    name : str or list
        The record name(s) to parse by.  For relaxed crystal records, the
        names should correspond to the key.
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
    return_df : bool, optional
        If True, then the corresponding pandas.Dataframe of metadata
        will also be returned.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    return self.get_records('relaxed_crystal', local=local, remote=remote, name=name,
                            key=key, method=method, standing=standing, family=family,
                            parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
                            potential_LAMMPS_key=potential_LAMMPS_key,
                            potential_id=potential_id, potential_key=potential_key,
                            symbols=symbols, natoms=natoms, natypes=natypes,
                            return_df=return_df, verbose=verbose)

def get_relaxed_crystal(self, local=None, remote=None, name=None,
                        key=None, method='dynamic', standing='good',
                        family=None, parent_key=None, 
                        potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                        potential_id=None, potential_key=None,
                        symbols=None, natoms=None, natypes=None, keyword=None,
                        prompt=True, verbose=False):
    """
    Retrieves exactly one matching relaxed crystal from the database.
    
    Parameters
    ----------
    local : bool, optional
        Indicates if the local location is to be searched.  Default value
        matches the value set when the database was initialized.
    remote : bool, optional
        Indicates if the remote location is to be searched.  Default value
        matches the value set when the database was initialized.
    name : str or list
        The record name(s) to parse by.  For relaxed crystal records, the
        names should correspond to the key.
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
    prompt : bool
        If prompt=True (default) then a screen input will ask for a selection
        if multiple matching potentials are found.  If prompt=False, then an
        error will be thrown if multiple matches are found.
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """
    def promptfxn(df):
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

    return self.get_record('relaxed_crystal', local=local, remote=remote, name=name,
                           key=key, method=method, standing=standing, family=family,
                           parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
                           potential_LAMMPS_key=potential_LAMMPS_key,
                           potential_id=potential_id, potential_key=potential_key,
                           symbols=symbols, natoms=natoms, natypes=natypes,
                           prompt=prompt, promptfxn=promptfxn, verbose=verbose)

def download_relaxed_crystals(self, name=None, key=None, method='dynamic',
                              standing='good', family=None, parent_key=None, 
                              potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                              potential_id=None, potential_key=None,
                              symbols=None, natoms=None, natypes=None, keyword=None,
                              overwrite=False, verbose=False):
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
    verbose : bool, optional
        If True, info messages will be printed during operations.  Default
        value is False.
    """

    self.download_records('relaxed_crystal', name=name,
                          key=key, method=method, standing=standing, family=family,
                          parent_key=parent_key, potential_LAMMPS_id=potential_LAMMPS_id,
                          potential_LAMMPS_key=potential_LAMMPS_key,
                          potential_id=potential_id, potential_key=potential_key,
                          symbols=symbols, natoms=natoms, natypes=natypes,
                          overwrite=overwrite, verbose=verbose)