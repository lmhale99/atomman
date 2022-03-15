# coding: utf-8

# Standard Python imports
from typing import Union

# https://github.com/usnistgov/potentials
import potentials

class Database(potentials.Database):
    """
    Child of potentials.Database extended for interacting with structure
    and defect reference records.
    """
    # Class imports
    from ._crystal_prototype import (get_crystal_prototypes, get_crystal_prototype,
                                     retrieve_crystal_prototype,
                                     download_crystal_prototypes)

    from ._relaxed_crystal import (get_relaxed_crystals, get_relaxed_crystal,
                                   retrieve_relaxed_crystal,
                                   download_relaxed_crystals)

    from ._reference_crystal import (get_reference_crystals, get_reference_crystal,
                                     retrieve_reference_crystal,
                                     download_reference_crystals, fetch_mp_crystal,
                                     fetch_oqmd_crystal, fetch_mp_crystals,
                                     fetch_reference_crystal)

    def download_all(self,
                     status: Union[str, list, None] = None,
                     downloadfiles: bool = True,
                     overwrite: bool = False,
                     verbose: bool = False):
        """
        Downloads potential and structure-related records from the remote
        location to the local location.

        Parameters
        ----------
        status : str, list or None, optional
            Only potential_LAMMPS records with the given status(es) will be
            downloaded.  Allowed values are 'active' , 'superseded', and 'retracted'.
            If None (default) is given, then all potentials will be downloaded.
        downloadfiles : bool, optional
            If True, the parameter files associated with the potential_LAMMPS
            record will also be downloaded.
        overwrite : bool, optional
            Flag indicating if any existing local records with names matching
            remote records are updated (True) or left unchanged (False).  Default
            value is False.
        verbose : bool, optional
            If True, info messages will be printed during operations.  Default
            value is False.
        """
        super().download_all(status=status, downloadfiles=downloadfiles,
                             overwrite=overwrite, verbose=verbose)

        self.download_crystal_prototypes(overwrite=overwrite, verbose=verbose)
        self.download_relaxed_crystals(overwrite=overwrite, verbose=verbose)
        #self.download_reference_crystals(overwrite=overwrite, verbose=verbose)