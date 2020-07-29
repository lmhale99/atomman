import potentials

from ...tools import aslist

class Database(potentials.Database):
    """
    Child of potentials.Database extended for interacting with structure
    and defect reference records.
    """
    # Class imports
    from ._crystal_prototype import (crystal_prototypes, crystal_prototypes_df,
                            load_crystal_prototypes, _no_load_crystal_prototypes,
                            get_crystal_prototypes, get_crystal_prototype,
                            download_crystal_prototypes)

    from ._relaxed_crystal import (relaxed_crystals, relaxed_crystals_df,
                            load_relaxed_crystals, _no_load_relaxed_crystals,
                            get_relaxed_crystals, get_relaxed_crystal,
                            download_relaxed_crystals)

    def __init__(self, host=None, username=None, password=None, certification=None,
                 localpath=None, verbose=False, local=True, remote=True, 
                 load=False, status='active'):

        # Extract extra loads options
        newoptions = ['crystal_prototype', 'relaxed_crystal']
        if load is True:
            newload = newoptions
        elif load is False:
            newload = []
        else:
            load = aslist(load)
            newload = []
            for key in newoptions:
                try:
                    newload.append(load.pop(load.index(key)))
                except:
                    pass
        
        # Call parent init
        super().__init__(host=host, username=username, password=password, certification=certification,
                 localpath=localpath, verbose=verbose, local=local, remote=remote, 
                 load=load, status=status)

        # Call new load options
        if 'crystal_prototype' in newload:
            self.load_crystal_prototypes(verbose=verbose)
        else:
            self._no_load_crystal_prototypes()

        # Call new load options
        if 'relaxed_crystal' in newload:
            self.load_relaxed_crystals(verbose=verbose)
        else:
            self._no_load_relaxed_crystals()

    def download_all(self, localpath=None, format='json', citeformat='bib',
                     indent=None, status='active', verbose=False,
                     get_files=True):
        """
        Downloads all records from the remote to localhost.

        Parameters
        ----------
        localpath : str, optional
            Path to a local directory where the records are to be copied to.
            If not given, will check localpath value set during object
            initialization.
        format : str, optional
            The file format to save the results locally as.  Allowed values are
            'xml' and 'json'.  Default value is 'json'.
        citeformat : str, optional
            The file format to save Citation records locally as.  Allowed
            values are 'xml', 'json', and 'bib'.  Default value is 'bib'.
        indent : int, optional
            The indentation spacing size to use for the locally saved record files.
            If not given, the JSON/XML content will be compact.
        verbose : bool, optional
            If True, info messages will be printed during operations.  Default
            value is False.
        status : str, list or None, optional
            Only potential_LAMMPS records with the given status(es) will be
            downloaded.  Allowed values are 'active' (default), 'superseded', and
            'retracted'.  If None is given, then all potentials will be downloaded.
        get_files : bool, optional
            If True, the parameter files associated with the potential_LAMMPS
            record will also be downloaded.
        
        Raises
        ------
        ValueError
            If no localpath, no potentials, invalid format, or records in a
            different format already exist in localpath.
        """
        self.download_citations(localpath=localpath, format=citeformat,
                                indent=indent, verbose=verbose)

        self.download_potentials(localpath=localpath, format=format,
                                 indent=indent, verbose=verbose)

        self.download_lammps_potentials(localpath=localpath, format=format,
                                        indent=indent, verbose=verbose,
                                        status=status, get_files=get_files)

        self.download_crystal_prototypes(localpath=localpath, format=format,
                                        indent=indent, verbose=verbose)