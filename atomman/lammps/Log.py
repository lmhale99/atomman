# coding: utf-8
# Standard Python libraries
import datetime
from typing import Optional, Union, List
import io

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# atomman imports
from ..tools import uber_open_rmode

class Simulation():
    """
    Object for representing the LAMMPS log output for a single MS/MD run
    """

    def __init__(self,
                 thermo: Optional[pd.DataFrame] = None,
                 performance: Optional[pd.DataFrame] = None):
        """
        Initializes a Simulation object with performance information.

        Parameters
        ----------
        thermo : pandas.DataFrame, optional
            The table of thermo data obtained from the simulation, if known.
        performance : pandas.DataFrame, optional
            The table of performance data obtained from the simulation, if known.
        """
        self.__thermo = None
        self.__performance = None
        self.__keys = []
        if thermo is not None:
            self.thermo = thermo
        if performance is not None:
            self.performance = performance

    def keys(self) -> list:
        """List of attribute keys that have been set"""
        return tuple(self.__keys)

    def __getitem__(self, key):
        """Accesses class attributes as dict keys""" 

        if key not in self.keys():
            raise KeyError(key)

        return getattr(self, key)

    def __iter__(self):
        """Iterates over set attribute keys"""
        for key in self.keys():
            yield key

    @property
    def thermo(self) -> Union[pd.DataFrame, None]:
        """pandas.DataFrame or None: The Simulation's thermo data"""
        return self.__thermo

    @property
    def performance(self) -> Union[pd.DataFrame, None]:
        """pandas.DataFrame or None: The Simulation's performance data"""
        return self.__performance

    @thermo.setter
    def thermo(self, value: pd.DataFrame):
        if isinstance(value, pd.DataFrame):
            self.__thermo = value
        else:
            self.__thermo = pd.DataFrame(value)
        if 'thermo' not in self.keys():
            self.__keys.append('thermo')

    @performance.setter
    def performance(self, value: pd.DataFrame):
        if isinstance(value, pd.DataFrame):
            self.__performance = value
        else:
            self.__performance = pd.DataFrame(value)
        if 'performance' not in self.keys():
            self.__keys.append('performance')

class Log():
    """Object for representing a LAMMPS log output"""

    def __init__(self,
                 log_info: Union[str, io.IOBase, None] = None):
        """
        Initializes a Log object.
        
        Parameters
        ----------
        log_info : str or file-like object, optional
            The LAMMPS log content to read in.  If None (default), then the
            Log object is created but empty.
        """

        # Initialize simulation properties
        self.__simulations = []
        self.__lammps_version = None
        self.__lammps_date = None

        # Read log data if supplied
        if log_info is not None:
            self.read(log_info)

    def read(self,
             log_info: Union[str, io.IOBase],
             append: bool = True):
        """
        Parses a LAMMPS screen output/log file.
        
        Parameters
        ----------
        log_info : str or file-like object
            The LAMMPS log content to read in.
        append : bool, optional
            Flag indicating if the content being read in is appended to the
            current data (True, default), or if it overwrites any existing
            saved content (False).
        """

        # Reset properties and values if append is False
        if append is False:
            self.__simulations = []
            self.__lammps_version = None
            self.__lammps_date = None

        # Strings found directly before and after run and minimize thermo data
        thermo_start_trigger = ['Memory usage per processor =',
                             'Per MPI rank memory allocation (min/avg/max) =']
        thermo_end_trigger = ['Loop time of']

        # Strings found directly before and after run and minimize performance data
        performance_start_trigger = ['MPI task timing breakdown']
        performance_start_trigger_old_version = ['Pair  time (%)']
        performance_end_trigger = ['Nlocal:']
        is_old_version = False

        # Handle file names, strings and open file-like objects equivalently
        with uber_open_rmode(log_info) as log_info:

            # Initialize parameters
            thermo_headers = []
            thermo_footers = []
            performance_headers = []
            performance_footers = []
            i = 0

            # For all lines in file/output
            for line in log_info:
                line = line.decode('UTF-8')

                # Skip blank lines
                if len(line.split()) == 0:
                    continue

                # Save the LAMMPS version information
                if line[:8] == 'LAMMPS (' and self.lammps_version is None:
                    self.__read_lammps_version(line)

                # Check for strings listed prior to run and minimize simulations
                if any([trigger in line for trigger in thermo_start_trigger]):
                    thermo_headers.append(i+1)

                # Check for strings listed after run and minimize simulations
                elif any([trigger in line for trigger in thermo_end_trigger]):
                    thermo_footers.append(i-1)

                # Check for strings listed prior to  performance data
                if any([trigger in line for trigger in performance_start_trigger]):
                    performance_headers.append(i+1)
                if any([trigger in line for trigger in performance_start_trigger_old_version]):
                    performance_headers.append(i)
                    is_old_version = True

                # Check for strings listed after performance data
                elif any([trigger in line for trigger in performance_end_trigger]):
                    performance_footers.append(i-1)

                i += 1

            # Add last line to footers for incomplete logs
            thermo_footers.append(i)

            # Reset file pointer
            log_info.seek(0)

            # Get number of Simulations already read
            j = len(self.simulations)

            # Read thermo data and create Simulations
            for header, footer in zip(thermo_headers, thermo_footers):
                self.__read_thermo(log_info, header, footer)

            # Read performance data for each Simulation
            for i in range(len(performance_footers)):
                header = performance_headers[i]
                footer = performance_footers[i]
                performance = self.__read_performance(log_info, header, footer, is_old_version)
                self.simulations[i+j].performance = performance

    def __read_lammps_version(self, line):
        """
        Subfunction for reading the LAMMPS version from the log file
        """
        month = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4,
                 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8,
                 'Sep': 9, 'Oct': 10,'Nov': 11,'Dec': 12}
        self.__lammps_version = line.strip()[8:-1]
        d = self.lammps_version.split('-')[0].split()
        self.__lammps_date = datetime.date(int(d[2]), month[d[1]], int(d[0]))

    def __read_thermo(self, log_info, header, footer):
        """
        Subfunction for reading the thermo data associated with a simulation run.
        """
        # Use pandas to read all thermo data at once
        thermo = pd.read_csv(log_info, header=header,
                                nrows=footer-header,
                                sep=r"\s+",
                                skip_blank_lines=True)

        # Reset file pointer
        log_info.seek(0)

        # Append simulation results
        self.__simulations.append(Simulation(thermo=thermo))

    def __read_performance(self, log_info, header, footer, is_old_version):
        """
        Subfunction for reading the performance data associated with a simulation run.
        """

        # Use pandas to read all performance data at once
        if not is_old_version:
            performance = pd.read_csv(log_info, header=header,
                                nrows=footer-header,
                                sep = '|',
                                skip_blank_lines=True)
            performance = performance.drop([0])
            performance.rename(columns=lambda x: x.strip(), inplace=True)
            performance = performance.set_index('Section')
            performance = performance.replace(r'^\s*$',0.0,regex=True)
            performance = performance.astype(float)

        else:
            performance = pd.read_csv(log_info, header=header,
                                nrows=footer-header,
                                sep = '=',
                                skip_blank_lines=True)
            performance = performance.columns.to_frame().T.append(performance,ignore_index=True)
            performance.rename(columns=lambda x: x.strip(), inplace=True)
            performance = performance.replace(')','')
            performance.columns = ['Section','time']
            performance = performance.set_index('Section')
            performance[['avg. Time','percentage']] = performance.time.str.split('(',expand=True)
            performance[['%','symbol']] = performance.percentage.str.split(')',expand=True)
            del performance['time']
            del performance['symbol']
            del performance['percentage']

        # Reset file pointer
        log_info.seek(0)

        # Return performance data
        return performance

    @property
    def simulations(self) -> List[dict]:
        """list of dict: parsed data for each separate LAMMPS run/minimize action"""
        return self.__simulations

    @property
    def lammps_version(self) -> str:
        """str : The LAMMPS version used."""
        return self.__lammps_version

    @property
    def lammps_date(self) -> datetime.date:
        """datetime.date : The date associated with the LAMMPS version used."""
        return self.__lammps_date

    def flatten(self,
                style: str = 'last',
                firstindex: Optional[int] = None,
                lastindex: Optional[int] = None) -> Simulation:
        """
        Combines multiple simulations into one.

        Parameters
        ----------
        style : str, optional
            Specifies how duplicate time step values are handled, i.e which
            values to keep:
            'first' uses the values from the earliest simulation.
            'last' uses the values from the latest simulation (default).
            'all' uses all reported lines including ones with duplicate time
            steps.
        firstindex : int or None, optional
            The leading list range index to limit which simulations are included
            in the merge, i.e. simulations[firstindex:lastindex].
        lastindex : int or None, optional
            The trailing list range index to limit which simulations are included
            in the merge, i.e. simulations[firstindex:lastindex].

        Returns
        -------
        Simulation
            A Simulation object containing the condensed thermo data.
        """
        simulations = self.simulations[firstindex:lastindex]

        # Check that all simulations with thermo data have step values
        for sim in simulations:
            if sim.thermo is not None and len(sim.thermo) > 0:
                assert 'Step' in sim.thermo, 'All simulation thermos must have Step key in order to flatten'

        # Combine the data into merged_df
        merged_df = simulations[0].thermo

        dtypes = {}
        for sim in simulations[1:]:
            thermo = sim.thermo

            if thermo is None:
                continue

            if style == 'first':
                merged_df = pd.concat([merged_df, thermo[thermo.Step > merged_df.Step.max()]], ignore_index=True)

            elif style == 'last':
                merged_df = pd.concat([merged_df[merged_df.Step < thermo.Step.min()], thermo], ignore_index=True)

                # Identify dtype of each column based on what it is in the last simulation where it appears
                for key in thermo.keys():
                    dtypes[key] = thermo[key].dtype

            elif style == 'all':
                merged_df = pd.concat([merged_df, thermo], ignore_index=True)

            else:
                raise ValueError('Unsupported style')

            # Fix for incomplete log values in crashed runs possibly altering data types
            for key in dtypes:
                try:
                    # Try setting dtype of column to match what it is in the final simulation run.
                    merged_df[key] = np.asarray(merged_df[key], dtype=dtypes[key])
                except ValueError:
                    pass

        return Simulation(thermo=merged_df)
