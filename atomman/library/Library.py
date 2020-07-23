# coding: utf-8
# Standard Python libraries
from pathlib import Path
import uuid

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/atomman
import atomman as am

import potentials

import requests

from .. import Settings
from ..tools import aslist

class Library():
    """
    Class for interacting with the iprPy library directory.
    """
    def __init__(self, directory=None, load=False, local=True, remote=True):
        """
        Class initializer

        Parameters
        ----------
        directory : str or Path, optional
            The path to the library directory to interact with.  If not given,
            will use the library directory that has been set for iprPy.
        """
        if directory is None:
            self.__directory = Settings().library_directory
        else:
            self.__directory = Path(directory).resolve()

        self.__potdb = potentials.Database(localpath=self.directory,
                                           local=local, remote=remote,
                                           load=load)

    @property
    def directory(self):
        return self.__directory

    @property
    def potdb(self):
        return self.__potdb

    @property
    def all_ref_styles(self):
        """list : all reference record styles hosted at potentials.nist.gov"""
        return ['crystal_prototype', 'dislocation', 'free_surface', 'point_defect', 'stacking_fault', 'potential_LAMMPS']

    def download_refs(self, style=None, status='active'):
        """
        Downloads reference records from potentials.nist.gov to the library.
        Note: this will overwrite any local copies of records with matching
        names.  If you made changes to library files, be sure to save them
        with a different name.

        Parameters
        ----------
        style : str or list, optional
            The reference style(s) to download.  If not given, all reference
            style will be downloaded.
        status : str, list or None, optional
            Only the potential_LAMMPS records with the given status(es) will
            be downloaded.  Allowed values are 'active' (default),
            'superseded', and 'retracted'.  If set to None, all hosted
            potential_LAMMPS will be downloaded.
        """          
        # Get all_ref_styles if none are specified
        if style is None:
            style = self.all_ref_styles
        style = aslist(style)
        
        for s in style:
            if s == 'potential_LAMMPS':
                self.potdb.download_lammps_potentials(format='json', indent=4, status=status, verbose=True)
            else:
                self.potdb.download_records(s, format='json', indent=4, verbose=True)
        
    def get_potentials(self, id=None, key=None, potid=None, potkey=None,
                    status='active', pair_style=None, element=None,
                    symbol=None, verbose=False, get_files=False):
        """
        Gets LAMMPS potentials from the iprPy library or by downloading from
        potentials.nist.gov if local copies are not found.
        
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
            
        Returns
        -------
        Potential.LAMMPSPotential
            The potential object to use.
        """
        if self.potdb.lammps_potentials is None:
            self.potdb.load_lammps_potentials()
            
        return self.potdb.get_lammps_potentials(id=id, key=key, potid=potid, potkey=potkey,
                                                status=status, pair_style=pair_style, element=element,
                                                symbol=symbol, verbose=verbose, get_files=get_files)

    def get_potential(self, id=None, key=None, potid=None, potkey=None,
                    status='active', pair_style=None, element=None,
                    symbol=None, verbose=False, get_files=False):
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
            
        Returns
        -------
        Potential.LAMMPSPotential
            The potential object to use.
        """
        
        return self.potdb.get_lammps_potential(id=id, key=key, potid=potid, potkey=potkey,
                                               status=status, pair_style=pair_style, element=element,
                                               symbol=symbol, verbose=verbose, get_files=get_files)