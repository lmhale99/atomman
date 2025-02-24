# coding: utf-8
# Standard Python libraries
import os
from typing import Optional, List, Tuple

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

from .Log import Log

class NEBLog(object):
    
    def __init__(self,
                 neblog: str = 'log.lammps',
                 replicalogs: str = 'log.lammps.*',
                 rootdir: Optional[str] = None):
        """
        Initializes reader for LAMMPS NEB calculation log files.
        
        Parameters
        ----------
        neblog : str, optional
            The filename for the primary LAMMPS NEB log file. Default value is
            'log.lammps'.
        replicalogs : str, optional
            The filename with wildcard for the LAMMPS NEB replica log files.
            Default value is 'log.lammps.*'
        rootdir : str, optional
            The root directory to use with respect to the neblog and
            replicalogs parameters.  Default value is None, which will assume
            the current working directory for relative paths.        
        """
        
        self.load(neblog=neblog, replicalogs=replicalogs,  rootdir=rootdir)

    @property
    def nreplicas(self) -> int:
        """int: The number of replicas"""
        return self.__nreplicas
    
    @property
    def minrun(self) -> pd.DataFrame:
        """pandas.DataFrame: The NEB log data for the minimization steps."""
        return self.__minrun
    
    @property
    def climbrun(self) -> pd.DataFrame:
        """pandas.DataFrame: The NEB log data for the barrier climb steps."""
        return self.__climbrun
    
    @property
    def logs(self) -> List[Log]:
        """list of atomman.lammps.log: The LAMMPS log files for each replica."""
        return self.__logs
    
    def load(self,
             neblog: str = 'log.lammps',
             replicalogs: str = 'log.lammps.*',
             rootdir: Optional[str] = None):
        """
        Loads LAMMPS NEB calculation log file data.
        
        Parameters
        ----------
        neblog : str, optional
            The filename for the primary LAMMPS NEB log file. Default value is
            'log.lammps'.
        replicalogs : str, optional
            The filename with wildcard for the LAMMPS NEB replica log files.
            Default value is 'log.lammps.*'
        rootdir : str, optional
            The root directory to use with respect to the neblog and
            replicalogs parameters.  Default value is None, which will assume
            the current working directory for relative paths.        
        """
        if rootdir is not None:
            neblog = os.path.join(rootdir, neblog)
            replicalogs = os.path.join(rootdir, replicalogs)
        replicalogs = replicalogs.replace('*', '%i')

        # First pass
        with open(neblog) as f:
            lines = f.readlines()

            # Count terms to identify number of replicas
            terms = lines[3].split()
            self.__nreplicas = int((len(terms) - 9) / 2)
            
            column_names = 'Step MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR RDT'.split()
            for i in range(self.nreplicas):
                column_names.append('RD%i' %(i+1))
                column_names.append('PE%i' %(i+1))

            climb_start = None
            nrows = None
            for i in range(3, len(lines)):
                if lines[i][:4] == 'Step':
                    climb_start = i + 1
                    nrows = i - 4

        # Second pass
        self.__minrun = pd.read_csv(neblog, names=column_names, skiprows=3, nrows=nrows, sep=r"\s+", skip_blank_lines=True)
        self.__climbrun = pd.read_csv(neblog, names=column_names, skiprows=climb_start, sep=r"\s+", skip_blank_lines=True)
        
        self.__logs = []
        for i in range(self.nreplicas):
            self.__logs.append(Log(replicalogs % (i) ))
       
    def get_neb_path(self, step: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Retrieves the reaction coordinates and corresponding potential
        energies for a given simulation step.
        
        Parameters
        ----------
        step : int
            The run step to retrieve values for.
            
        Returns
        -------
        reaction_coordinates : numpy.ndarray
            The reaction coordinates
        potential_energies : numpy.ndarray
            The potential energies
        """
        reaction_coordinate = []
        potential_energy = []
        joined = pd.concat([self.minrun, self.climbrun])
        match = joined[joined.Step == step]
        if len(match) > 0:
            match = match.iloc[0]
        else:
            raise ValueError('Step value not found in NEB log file')
        key = 'RD%i'
        for i in range(1, self.nreplicas+1):
            reaction_coordinate.append(match[key%i])    

        for replica in range(self.nreplicas):
            match = None    
            for sim in self.logs[replica].simulations:
                try:
                    match = sim['thermo'][sim['thermo'].Step == step].iloc[0]
                except:
                    continue
                else:
                    break
            if match is None:
                raise ValueError('Step value not found in log.lammps.%i file' %replica)
            potential_energy.append(match.PotEng)

        reaction_coordinates = np.array(reaction_coordinate)
        potential_energies = np.array(potential_energy) - potential_energy[0]

        return reaction_coordinates, potential_energies

    def get_barrier(self, reverse: bool = False) -> float:
        """
        Returns the barrier energy calculated from the final NEB simulation step.
        
        Parameters
        ----------
        reverse : bool, optional
            Indicates if the energy barrier returned is the forward barrier
            relative to the first replica (False, default), or is the reverse
            barrier relative to the last replica (True).
            
        Returns
        -------
        float
            The energy barrier relative to one of the endpoint replicas.
        """
        potential_energies = self.get_neb_path(self.climbrun.Step.values[-1])[1]
        
        if reverse is False:
            reference = potential_energies[0]
        elif reverse is True:
            reference = potential_energies[-1]
        else:
            raise ValueError('Parameter reverse must be bool')
        
        return np.max(potential_energies) - reference