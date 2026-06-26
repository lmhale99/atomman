from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import shlex
import shutil

from typing import Optional, Union
import datetime

import pandas as pd

# Local imports
import atomman.unitconvert as uc
from atomman.typing import lammps, lammpspotential, unitfloat
from . import checkversion, style, run, restart_check, read_logs

if TYPE_CHECKING:
    from .. import System

class dummy_command_wrapper:
    """Wrapper class to enable using 'lmp.xxx("args")' instead of 'lmp.command("xxx args")'"""
    def __init__(self, lmp):
        self.lmp = lmp
        self.auto_flush = False

    def __dir__(self):
        return sorted(set(['angle_coeff', 'angle_style', 'atom_modify', 'atom_style', 'atom_style',
        'bond_coeff', 'bond_style', 'boundary', 'change_box', 'communicate', 'compute',
        'create_atoms', 'create_box', 'delete_atoms', 'delete_bonds', 'dielectric',
        'dihedral_coeff', 'dihedral_style', 'dimension', 'dump', 'fix', 'fix_modify',
        'group', 'improper_coeff', 'improper_style', 'include', 'kspace_modify',
        'kspace_style', 'lattice', 'mass', 'minimize', 'min_style', 'neighbor',
        'neigh_modify', 'newton', 'nthreads', 'pair_coeff', 'pair_modify',
        'pair_style', 'processors', 'read', 'read_data', 'read_restart', 'region',
        'replicate', 'reset_timestep', 'restart', 'run', 'run_style', 'thermo',
        'thermo_modify', 'thermo_style', 'timestep', 'undump', 'unfix', 'units',
        'variable', 'velocity', 'write_restart'] + self.lmp.available_styles("command")))

    def __getattr__(self, name):
        """
        This method is where the Python 'magic' happens. If a method is not
        defined by the class command_wrapper, it assumes it is a LAMMPS command. It takes
        all the arguments, concatinates them to a single string, and executes it using
        :py:meth:`lammps.command`.

        LAMMPS commands that accept callback functions (such as fix python/invoke)
        can be passed functions and lambdas directly. The first argument of such
        callbacks will be an lammps object constructed from the passed LAMMPS
        pointer.

        :return: line or list of lines of output, None if no output
        :rtype: list or string
        """
        def handler(*args, **kwargs):
            cmd_args = [name] + [str(x) for x in args]

            # Python 3.6+ maintains ordering of kwarg keys
            for kw, arg in kwargs.items():
                cmd_args.append(kw)
                if isinstance(arg, bool):
                    cmd_args.append("true" if arg else "false")
                else:
                    cmd_args.append(str(arg))

            cmd = ' '.join(cmd_args)
            self.lmp.command(cmd)
            if self.auto_flush:
                self.lmp.flush_buffers()
        return handler


class DummyLAMMPS():
    """Dummy LAMMPS interfafce that essentially does nothing when any attributes are accessed"""

    def __init__(self,name='',cmdargs=None,ptr=None,comm=None):
        self._cmd = dummy_command_wrapper(self)

    @property
    def callables(self) -> list:
        """List of all callable attributes of lammps.lammps"""
        return ['addstep_compute', 'addstep_compute_all', 'available_ids', 'available_plugins', 'available_styles', 'c_bigint',
                'c_imageint', 'c_tagint', 'clearstep_compute', 'close', 'command', 'commands_list', 'commands_string',
                'create_atoms', 'create_molecule', 'decode_image_flags', 'encode_image_flags', 'error', 'eval', 'expand',
                'extract_atom', 'extract_atom_datatype', 'extract_atom_size', 'extract_box', 'extract_compute', 'extract_fix',
                'extract_global', 'extract_global_datatype', 'extract_pair', 'extract_pair_dimension', 'extract_setting',
                'extract_variable', 'file', 'finalize', 'find_compute_neighlist', 'find_fix_neighlist', 'find_pair_neighlist',
                'fix_external_get_force', 'fix_external_set_energy_global', 'fix_external_set_energy_peratom', 'fix_external_set_vector',
                'fix_external_set_vector_length', 'fix_external_set_virial_global', 'fix_external_set_virial_peratom',
                'flush_buffers', 'force_timeout', 'gather', 'gather_angles', 'gather_atoms', 'gather_atoms_concat',
                'gather_atoms_subset', 'gather_bonds', 'gather_concat', 'gather_dihedrals','gather_impropers', 'gather_subset',
                'get_gpu_device_info', 'get_mpi_comm', 'get_natoms', 'get_neighlist', 'get_neighlist_element_neighbors',
                'get_neighlist_size', 'get_os_info', 'get_thermo', 'has_id', 'has_package', 'has_style', 'last_thermo',
                'map_atom', 'reset_box', 'scatter', 'scatter_atoms', 'scatter_atoms_subset', 'scatter_subset',
                'set_fix_external_callback', 'set_internal_variable', 'set_show_error', 'set_string_variable',
                'set_variable', 'version']
    
    @property
    def properties(self) -> list:
        """List of all non-callable, i.e. properties of lammps.lammps"""
        return ['accelerator_config', 'callback', 'cmd', 'comm', 'has_curl_support', 'has_exceptions', 'has_ffmpeg_support',
                'has_gpu_device', 'has_gzip_support', 'has_jpeg_support', 'has_mpi4py', 'has_mpi_support', 'has_png_support',
                'installed_packages', 'ipython', 'is_running', 'last_thermo_step', 'lib', 'lmp', 'numpy', 'opened']

    def __getattr__(self, name):
        """
        dummy switch for functions or properties
        """
        if name in self.callables:
            return self.dummyfxn
        elif name in self.properties:
            return None
        else:
            raise AttributeError(f"'lammps' object has no attribute '{name}'")
    
    def dummyfxn(self, *args, **kwargs):
        return None
    
    def available_styles(self, *args, **kwargs):
       return []
    
    @property
    def cmd(self):
        return self._cmd


class ExtendLAMMPS():
    """Common attributes used by LAMMPSLIB and LAMMPSEXE objects"""

    def __init__(self,
                 potential: Optional[lammpspotential] = None,
                 logfile: str = 'log.lammps'):
        
        # Set islib and versiondate for a LAMMPSLIB object
        if isinstance(self, lammps.lammps):
            self.__islib = True

        # Set islib and versiondate for a LAMMPSEXE object
        elif hasattr(self, 'lammps_command'):
            self.__islib = False

        else:
            raise ValueError('no lammps interface or lammps_command found!')
        
        # Get version information
        version_info = checkversion(self)
        self.__versionstr = version_info['version']
        self.__versiondate = version_info['date']

        # Set potential
        self.potential = potential
        
        # Initialize script
        self.script = ''

        # Set initial log details
        self.__logfile = logfile

        # Set default lognumber (indicates restart status)
        self.__lognumber = 0

    # ----------------------New Properties ---------------------------------- #

    @property
    def islib(self) -> bool:
        """bool: True if is LAMMPSLIB, False if is LAMMPSEXE"""
        return self.__islib
    
    @property
    def lognumber(self) -> int:
        """int: Number of log files from previous restart runs to read"""
        return self.__lognumber

    @property
    def script(self) -> str:
        """str: The current LAMMPS script of commands"""
        return self.__script
    
    @script.setter
    def script(self, val):
        assert isinstance(val, str)
        self.__script = val

    @property
    def versiondate(self) -> datetime.date:
        """datetime.date : The LAMMPS version date"""
        return self.__versiondate
    
    @property
    def versionstr(self) -> str:
        """str : The LAMMPS version string"""
        return self.__versionstr
    
    @property
    def potential(self) -> lammpspotential:
        """PotentialLAMMPS or PotentialLAMMPSKIM : The potential object"""
        if self.__potential is None:
            raise ValueError('No potential set!')
        return self.__potential
    
    @potential.setter
    def potential(self, val: Optional[lammpspotential]):
        if val is None:
            self.__potential = None
            self.__unitsdict = None
        elif isinstance(val, (lammpspotential)):
            self.__potential = val
            self.__unitsdict = style.unit(val.units)
        else:
            raise TypeError('potential must be a PotentialLAMMPS or PotentialLAMMPSKIM object, or None')

    @property
    def unitsdict(self) -> dict:
        """dict : The units of values associated with the LAMMPS units setting"""
        if self.__unitsdict is None:
            raise ValueError('No potential set!')
        return self.__unitsdict
    
    @property
    def logfile(self) -> str:
        """str : The name/path of the LAMMPS log file. 'none' means no log file"""
        return self.__logfile

    # ---------------------Commands to script-------------------------------- #

    def command(self, cmd: str):
        """copy command string to script"""
        self.script += cmd + '\n'

    def commands_string(self, multicmd: str):
        """Copy multi-line command string to script"""
        self.script += multicmd + '\n'

    def commands_list(self, cmdlist: list):
        """Copy list of command strings to script"""
        self.script += '\n'.join(cmdlist) + '\n'

    def file(self, path: str):
        """Copy file contents to script"""
        with open(path) as f:
            self.script += f.read() + '\n'

    # ---------------------------General operations-------------------------- #

    def new_system_no_atoms(self,
                            system: System,
                            potential: Optional[lammpspotential] = None,
                            tilt_large: bool = False,
                            logfile: Optional[str] = None):
        """
        Generates the box and pair commands for setting up a new LAMMPS
        calculation but DOES NOT DEFINE ANY ATOMS!.  This allows for similar
        setup as the other "new_system" calls but leaves atom creation to
        be defined elsewhere.

        Parameters
        ----------
        system : atomman.System
            A system collecting box, pbc and symbols info to pass to LAMMPS.
            NOTE ANY INCLUDED ATOMS WILL BE IGNORED!
        potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
            Potential information to pass in.  Optional here if set
            during LAMMPS object initialization.
        tilt_large : bool, optional
            For older LAMMPS, will add the "box tilt large" command to
        logfile : str or None, optional
            Specifies a name/path for where LAMMPS log information is to be
            saved.  Default value of None will not change previous log file
            settings.
        """
        if potential is None:
            potential = self.potential
        else:
            self.potential = potential

        if logfile is not None and logfile != self.logfile:
            self.cmd.log(logfile)
            self.__logfile = logfile

        # Dump system and potential info to lmp object without atoms
        system.dump('lammps_lib', lmp=self, potential=self.potential,
                    lammps_date=self.versiondate, tilt_large=tilt_large, no_atoms=True)

    def new_system_from_data_file(self,
                                  system: System,
                                  filename: str = 'init.dat',
                                  potential: Optional[lammpspotential] = None,
                                  tilt_large: bool = False,
                                  usefiles: bool = False,
                                  logfile: Optional[str] = None,
                                  clear: bool = True):
        """
        (Re)sets the box, atoms and pair commands for setting up a new LAMMPS
        calculation based on a system object, which is passed to LAMMPS either
        directly or by saving it to a LAMMPS data file which LAMMPS then reads
        in.

        Parameters
        ----------
        system : atomman.System
            The system information to pass to LAMMPS
        filename : str, optional
            The name/path of the LAMMPS data file where system information
            is saved to and loaded in by LAMMPS.  This is always created
            by LAMMPSEXE objects and is optionally created by LAMMPSLIB
            objects based on the usefiles setting.
        potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
            Potential information to pass in.  Optional here if set
            during LAMMPS object initialization.
        tilt_large : bool, optional
            For older LAMMPS, will add the "box tilt large" command to
            keep LAMMPS from throwing an error.
        usefiles : bool, optional
            Setting this to True forces a LAMMPSLIB object to create and
            load the data file rather than passing in the system information
            directly.
        logfile : str or None, optional
            Specifies a name/path for where LAMMPS log information is to be
            saved.  Default value of None will not change previous log file
            settings.
        clear : bool, optional
            If True (default), will do a "clear" command in LAMMPS before
            defining the system.  Setting this to False allows for specifying
            settings that need to come prior to the system definition.
        """
        if potential is None:
            potential = self.potential
        else:
            self.potential = potential

        if logfile is not None and logfile != self.logfile:
            self.cmd.log(logfile)
            self.__logfile = logfile
        
        if self.islib and not usefiles:
            # pass system and pair info directly to LAMMPS
            system.dump('lammps_lib', lmp=self, potential=self.potential,
                        lammps_date=self.versiondate, tilt_large=tilt_large,
                        clear=clear)
        else:
            # Build init.dat, and create load and pair commands
            system_info = system.dump('atom_data', f=filename,
                                      potential=self.potential,
                                      tilt_large=tilt_large,
                                      lammps_date=self.versiondate,
                                      triclinic=True)
            if clear:
                self.cmd.clear()
            self.commands_string(system_info)

    def new_system_from_restart(self,
                                system: System,
                                filename: str,
                                potential: Optional[lammpspotential] = None,
                                tilt_large: bool = False,
                                usefiles: bool = False,
                                logfile: Optional[str] = None,
                                clear: bool = True):
        """
        (Re)sets the box, atoms and pair commands for setting up a new LAMMPS
        calculation based on a restart file or an equivalent system object.

        This supports two primary uses:
        1. If restarting incomplete simulations, set usefiles=True to
           always read from the restart file.
        2. If reverting to an earlier simulation state, usefiles=False
           allows for LAMMPSLIB objects to update the system info directly
           rather than reading in from the restart file.

        Parameters
        ----------
        system : atomman.System
            The system information to pass to LAMMPS.  If usefiles=True,
            only symbols information will be extracted from the system.
            If usefiles=False, the system info will be passed to a
            LAMMPSLIB object directly instead of usig the restart file.
        filename : str
            The restart filename value to use with the LAMMPS read_restart
            command.  This can include wildcards.
        potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
            Potential information to pass in.  Optional here if set
            during LAMMPS object initialization.
        tilt_large : bool, optional
            For older LAMMPS, will add the "box tilt large" command to
            keep LAMMPS from throwing an error.
        usefiles : bool, optional
            Setting this to True forces a LAMMPSLIB object to use the
            restart file rather than passing in the system information
            directly.
        logfile : str or None, optional
            Specifies a name/path for where LAMMPS log information is to be
            saved.  Default value of None will not change previous log file
            settings.
        """
        if potential is None:
            potential = self.potential
        else:
            self.potential = potential

        if logfile is not None and logfile != self.logfile:
            self.cmd.log(logfile)
            self.__logfile = logfile
        
        if self.islib and not usefiles and system is not None:
            # pass system and pair info directly to LAMMPS
            system.dump('lammps_lib', lmp=self, potential=self.potential,
                        tilt_large=tilt_large)
        else:
            # Create load restart and pair commands
            system_info = potential.pair_restart_info(filename = filename,
                                                      symbols = system.symbols,
                                                      lammps_date = self.versiondate,
                                                      tilt_large = tilt_large)
            if clear:
                self.cmd.clear()
            self.commands_string(system_info)

    def loop(self,
             name: str,
             number: int,
             usefiles: bool = False,
             ):
        """
        Common handling of loops for both exe and lib variations of LAMMPS.
        
        For exe, this adds the LAMMPS variable loop, label, next, and jump
        commands around the embedded commands.
        
        For lib, this updates the LAMMPS variable value and yields it allowing
        looping to be handled in Python.

        NOTE: The looped values are consistent with LAMMPS variable loop and
        range from 1 to number!

        Parameters
        ----------
        name : str
            A name for the looping variable.
        number : int
            The max iteration number.
        usefiles : bool
            If True, will always use LAMMPS looping even with LAMMPS libs.

        Yields
        ------
        name : str
            The variable name is yielded if exe or usefiles=True
        n : int
            The iteration value (1 to number) if lib and usefiles=False
        """
        # Use LAMMPS looping
        if not self.islib or usefiles:
            # Open LAMMPS loop
            self.cmd.variable(name, 'loop', number)
            self.cmd.label(f'loop_{name}')

            yield name

            # Close LAMMPS loop
            self.commands_string(f'\n# Close loop {name}')
            self.cmd.next(name)
            self.cmd.jump('SELF', f'loop_{name}')

        # Use Python looping
        else:
            for n in range(1, number+1):

                # Update variable value in LAMMPS and Python
                self.cmd.variable(name, 'equal', n)
                yield n

    def restart_check(self,
                      logfile: str = 'log.lammps',
                      restart_filename: Optional[str] = None,
                      ) -> bool:
        """
        Checks for existing restart and log files indicating that incomplete
        calculations should be continued.  Restart state can be determined
        by this method's output or by checking the value of lognumber after
        calling.

        Parameters
        ----------
        logfile : str, optional
            The LAMMPS log file name to seach for.  If found, it will be
            renamed to avoid being written over.  Default value is
            "log.lammps".
        restart_filename : str
            The restart file name to search for.  Can include glob-recognized
            wildcards.
        

        Returns
        -------
        bool
            True if restart files exist, False otherwise.
        """
        logname, logext, lognum = restart_check(logfile, restart_filename)

        self.__lognumber = lognum

        return lognum > 0

    def end_and_get_log(self,
                        scriptname: Optional[str] = None,
                        partition: Optional[str] = None,
                        return_log: bool = True):
        """
        Clears the saved script, optionally saving it.

        Parameters
        ----------
        scriptname : str or None, optional
            A file name where the compiled script will be saved.  Default
            value of None will not save it.
        
        Returns
        -------
        atomman.lammps.Log
            The interpreted log/screen output.
        """
        # Save the script 
        if scriptname is not None:
            with open(scriptname, 'w') as f:
                f.write(self.script)
            script = None
            script_name = scriptname
        else:
            script = self.script
            script_name = None
        
        if not self.islib:

            # LAMMPSEXE runs the script using a LAMMPS executable
            screen = self.logfile == 'none'
            log = run(self.lammps_command, script=script, script_name=script_name,
                      mpi_command=self.mpi_command, logfile=self.logfile,
                      restart_lognum = self.lognumber, partition=partition,
                      return_log=return_log, screen=screen)
        
        elif self.logfile != 'none' and return_log:
            # LAMMPSLIB reads in log files if they exist
            log = read_logs(self.logfile, lognum=self.lognumber)
            
        else:
            # If no log file, expect that matching contents are read in elsewhere...
            log = None

        # Clear and reset
        self.script = ''
        self.__lognumber = 0
        
        return log

    def set_thermo_units(self, thermo: Union[dict, pd.Series, pd.DataFrame]):
        """
        Converts all standard thermo terms from LAMMPS units into working units

        Parameters
        ----------
        thermo : dict, pandas.Series, pandas.DataFrame
            The thermo output data.
        
        Returns
        -------
        unconverted : list
            The keys in thermo that were not recognized as standard thermo keys
            and therefore were not converted.
        """
        # Define conversion units based on thermo keys and LAMMPS units setting
        thermo_units = {}
        
        timevals = ['Dt', 'Time']
        for timeval in timevals:
            thermo_units[timeval] = self.unitsdict['time']

        energyvals = ['PotEng', 'KinEng', 'TotEng', 'E_vdwl', 'E_coul', 'E_pair',
                      'E_bond', 'E_angle', 'E_dihed', 'E_impro', 'E_mol', 'E_long',
                      'E_tail', 'Enthalpy', 'Ecouple', 'Econserve']
        for energyval in energyvals:
            thermo_units[energyval] = self.unitsdict['energy']
            
        pressurevals = ['Press', 'Pxx', 'Pyy', 'Pzz', 'Pxy', 'Pxz', 'Pyz']
        for pressureval in pressurevals:
            thermo_units[pressureval] = self.unitsdict['pressure']

        lengthvals = ['Xlo',  'Xhi', 'Ylo', 'Yhi', 'Zlo', 'Zhi', 'Xy', 'Xz', 'Yz',
                      'Avecx', 'Avecy', 'Avecz', 'Bvecx', 'Bvecy', 'Bvecz',
                      'Cvecx', 'Cvecy', 'Cvecz', 'Lx', 'Ly', 'Lz',
                      'Cella' 'Cellb', 'Cellc']
        for lengthval in lengthvals:
            thermo_units[lengthval] = self.unitsdict['length']

        forcevals = ['Fmax', 'Fnorm']
        for forceval in forcevals:
            thermo_units[forceval] = self.unitsdict['force']

        # Unique conversion units
        thermo_units['Volume'] = self.unitsdict['length']+'^3'
        thermo_units['Density'] = self.unitsdict['density']

        # List of known thermo keys that conversions are not done for
        no_convert = ['Step', 'Elapsed', 'Elaplong', 'CPU', 'T/CPU', 'S/CPU', '%CPU',
                      'CPULeft', 'Part', 'TimeoutLeft', 'Atoms', 'Temp',
                      'Xlat', 'Ylat', 'Zlat', 'CellAlpha', 'CellBeta', 'CellGamma',
                      'Bonds', 'Angles', 'Diheds', 'Impros', 'Nbuild', 'Ndanger']
    
        unconverted = []
        for thermo_key in thermo:
            if thermo_key in thermo_units:
                thermo[thermo_key] = uc.set_in_units(thermo[thermo_key], thermo_units[thermo_key])
            elif thermo_key not in no_convert:
                unconverted.append(thermo_key)

        return unconverted

    # -----------------------LAMMPS calculation chunks----------------------- #

    def define_thermo(self,
                      thermosteps: int = 0,
                      thermokeys: Union[str, list, None] = None):
        """
        Defines the LAMMPS script lines associated with setting up the thermo
        output

        Parameters
        ----------
        thermosteps : int, optional
            How often thermo data is to be generated.  Default value of 0
            only outputs starting and ending values for a simulation run.
        thermokeys : str, list or None, optional
            What thermo values to output, given as either a space-delimited
            string or a list.  Default value of None outputs step, potential
            energy and pressures.
        """
        # Handle thermokeys values
        if thermokeys is None:
            thermokeys = 'step pe pxx pyy pzz pyz pxz pxy'
        if isinstance(thermokeys, str):
            thermokeys = thermokeys.split(' ')

        # Build thermo command lines
        self.cmd.thermo(thermosteps)
        self.cmd.thermo_style('custom', *thermokeys)
        self.cmd.thermo_modify('format', 'float', '%.17e')

    def define_dump(self,
                    dumpsteps: int,
                    dumpkeys: Union[str, list, None] = None,
                    dumpfile: str = '*.dump'):
        """
        Defines the LAMMPS script lines associated with setting up dump
        output

        Parameters
        ----------
        dumpsteps : int, optional
            How often dump files are to be generated.
        dumpkeys : str, list or None, optional
            What atomic values to output, given as either a space-delimited
            string or a list.  Default value of None outputs id, type, x, y and z.
        dumpfile : str, optional
            The dump file name(s) to use.  Default value is '*.dump'.
        """
        # Handle dumpkeys values
        if dumpkeys is None:
            dumpkeys = 'id type x y z'
        if isinstance(dumpkeys, str):
            dumpkeys = dumpkeys.split(' ')

        self.cmd.dump('dumpit', 'all', 'custom', dumpsteps, dumpfile, *dumpkeys)
        self.cmd.dump_modify('dumpit', 'format', 'float', '%.17e')


    def define_restart(self,
                       restartsteps: int,
                       restartfile: str = '*.restart'):
        """
        Defines the LAMMPS script lines associated with setting up restart
        output

        Parameters
        ----------
        restartsteps : int, optional
            How often restart files are to be generated.
        restartfile : str, optional
            The restart file name(s) to use.  Default value is '*.restart'.
        """
        self.cmd.restart(restartsteps, restartfile)

    def define_nve(self,
                   timestep: Optional[float] = None):
        
        if timestep is None:
            timestep = style.timestep(self.potential.units)

        self.cmd.timestep(timestep)
        self.cmd.fix('nve', 'all', 'nve')

    def minimize(self,
                 etol: float = 0.0,
                 ftol: unitfloat = 0.0,
                 maxiter: int = 10000,
                 maxeval: int = 100000,
                 dmax: unitfloat = '0.01 angstrom'):

        # Convert values given with units if needed
        dmax = uc.set_in_units(dmax)
        ftol = uc.set_in_units(ftol)

        # Convert values to lammps units
        lammps_units = self.unitsdict
        dmax = uc.get_in_units(dmax, lammps_units['length'])
        ftol = uc.get_in_units(ftol, lammps_units['force'])

        # Define LAMMPS commands
        self.cmd.min_modify('dmax', dmax)
        self.cmd.minimize(etol, ftol, maxiter, maxeval)
    
    def run(self,
            runsteps: int,
            upto: bool = False):

        if upto:
            self.cmd.run(runsteps, 'upto')
        else:
            self.cmd.run()




class LAMMPSLIB(lammps.lammps, ExtendLAMMPS):

    def __init__(self,
                 name: str = '',
                 cmdargs: Optional[list] = None,
                 ptr = None,
                 comm = None,
                 potential: Optional[lammpspotential] = None):
        """
        Create a LAMMPSLIB object

        Parameters
        ----------
        name : str, optional
            The "machine" name to use.
            "mpi" loads 'liblammps_mpi.so, "" (default) loads liblammps.so.
        cmdargs : list or None, optional
            The list of any command line
            arguments to be passed to the :cpp:func:`lammps_open` function.
            The executable name is automatically added.
        ptr : pointer or None, optional
            A pointer to a LAMMPS C++ class
            instance when called from an embedded Python interpreter.  None
            (default) means load symbols from shared library.
        comm: mpi4py.MPI_Comm
            A communicator object for running
            LAMMPS in parallel. None (default) means use MPI_COMM_WORLD implicitly.
        potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
            A potential object to associate with the simulation.  Can be either
            given here or set later.
        """
        if cmdargs is not None:
            if '-l' in cmdargs:
                logfile = cmdargs[cmdargs.index('-l') + 1]
            elif '-log' in cmdargs:
                logfile = cmdargs[cmdargs.index('-log') + 1]
            else:
                logfile = 'log.lammps'
        else:
            logfile = 'log.lammps'
        
        # Call init of parent classes
        lammps.lammps.__init__(self, name, cmdargs, ptr, comm)
        ExtendLAMMPS.__init__(self, potential=potential, logfile=logfile)


    def command(self, cmd: str):
        """
        Process a single LAMMPS input command from a string.

        This is a wrapper around the :cpp:func:`lammps_command`
        function of the C-library interface.

        :param cmd: a single lammps command
        :type cmd:  string
        """
        # Call parent class functions
        ExtendLAMMPS.command(self,cmd)
        lammps.lammps.command(self,cmd)

    def commands_string(self, multicmd: str):
        """
        Process a block of LAMMPS input commands from a string.

        This is a wrapper around the
        :cpp:func:`lammps_commands_string`
        function of the C-library interface.

        :param multicmd: text block of lammps commands
        :type multicmd:  string
        """
        # Call parent class functions
        ExtendLAMMPS.commands_string(self,multicmd)
        lammps.lammps.commands_string(self,multicmd)

    def commands_list(self,cmdlist: list):
        """
        Process multiple LAMMPS input commands from a list of strings.

        This is a wrapper around the
        :cpp:func:`lammps_commands_list` function of
        the C-library interface.

        :param cmdlist: a single lammps command
        :type cmdlist:  list of strings
        """
        # Call parent class functions
        ExtendLAMMPS.commands_list(self,cmdlist)
        lammps.lammps.commands_list(self,cmdlist)

    def file(self, path: str):
        """
        Read LAMMPS commands from a file.

        This is a wrapper around the :cpp:func:`lammps_file` function of the C-library interface.
        It will open the file with the name/path `file` and process the LAMMPS commands line by line until
        the end. The function will return when the end of the file is reached.

        :param path: Name of the file/path with LAMMPS commands
        :type path:  string
        """
        # Call parent class functions
        ExtendLAMMPS.file(self,path)
        lammps.lammps.file(self,path)


class LAMMPSEXE(DummyLAMMPS, ExtendLAMMPS):

    def __init__(self,
                 lammps_command: str,
                 mpi_command: Optional[str] = None,
                 potential: Optional[lammpspotential] = None):
        """
        Creates a mock LAMMPS library interface that simply generates the
        associated LAMMPS command script to be executed later.

        Parameters
        ----------
        lammps_command : str or None, optional
            This is the executable to use plus optional command line arguments.
        mpi_command : str or None, optional
            The MPI executable and inline arguments to use when calling the
            LAMMPS script.  If None (default) then LAMMPS will run serially.
        potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
            A potential object to associate with the simulation.  Can be either
            given here or set later.
        """
        self.__lammps_command = lammps_command
        self.__mpi_command = mpi_command
        
        # Call inits of parent classes
        ExtendLAMMPS.__init__(self, potential=potential, logfile='log.lammps')
        DummyLAMMPS.__init__(self)

    @property
    def lammps_command(self) -> str:
        return self.__lammps_command
    
    @property
    def mpi_command(self) -> Optional[str]:
        return self.__mpi_command

    @mpi_command.setter
    def mpi_command(self, val: Optional[str]):
        if val is None:
            self.__mpi_command = None
        else:
            self.__mpi_command = str(val)

# Define extra typing hints
LAMMPSobj = Union[LAMMPSLIB, LAMMPSEXE] # atomman LAMMPS objects

def LAMMPS(*args: Union[str, LAMMPSobj],
           name: str = '',
           cmdargs: Optional[list] = None,
           ptr = None,
           comm = None,
           lammps_command: Optional[str] = None,
           mpi_command: Optional[str] = None,
           potential: Optional[lammpspotential] = None) -> LAMMPSobj:
    """
    Wrapper function to build either a LAMMPSEXE or LAMMPSLIB object based on
    the given inputs.

    Parameters
    ----------
    *args : str, LAMMPSEXE, or LAMMPSLIB, optional
        One single non-keyword argument is allowed, which if it is a string
        is interpreted as either a value for the name or lammps_command arguments
        depending on if it contains a path to an executable or not.  LAMMPSLIB
        and LAMMPSEXE objects are passed through.
    name : str, optional
        For LAMMPS library interfaces, this specifies the "machine" name to use.
        "mpi" loads 'liblammps_mpi.so, "" (default) loads liblammps.so.
    cmdargs : list or None, optional
        For LAMMPS library interfaces, this is the list of any command line
        arguments to be passed to the :cpp:func:`lammps_open` function.
        The executable name is automatically added.
    ptr : pointer or None, optional
        For LAMMPS library interfaces, this is a pointer to a LAMMPS C++ class
        instance when called from an embedded Python interpreter.  None
        (default) means load symbols from shared library.
    comm: mpi4py.MPI_Comm or None, optional
        For LAMMPS library interfaces, this is a communicator object for running
        LAMMPS in parallel. None (default) means use MPI_COMM_WORLD implicitly.
    lammps_command : str or None, optional
        For LAMMPS executables, this is the executable to use plus optional
        command line arguments.  If given, then a LAMMPSEXE object is returned
        and therefore cannot be combined with the above parameters.
    mpi_command : str or None, optional
        The MPI executable and inline arguments to use when calling the
        LAMMPS script.  If None (default) then LAMMPS will run serially.
    potential : PotentialLAMMPS, PotentialLAMMPSKIM or None, optional
        A potential object to associate with the simulation.  Can be either
        given here or set later.

    Returns
    -------
    LAMMPSEXE or LAMMPSLIB
        The LAMMPS interface object.
    """
    # Interpret args option
    if len(args) == 1:
        
        # Pass LAMMPSEXE through
        if isinstance(args[0], (LAMMPSEXE, LAMMPSLIB)):
            if (name != '' or cmdargs is not None or ptr is not None or comm is not None or
                lammps_command is not None):
                raise ValueError('LAMMPSEXE object already exists: cannot set parameters!')
            if mpi_command is not None:
                args[0].mpi_command = mpi_command
            if potential is not None:
                args[0].potential = potential
            return args[0]
        
       # Pass LAMMPSLIB through
        if isinstance(args[0], (LAMMPSEXE, LAMMPSLIB)):
            if (name != '' or cmdargs is not None or ptr is not None or comm is not None or
                lammps_command is not None or mpi_command is not None):
                raise ValueError('LAMMPSLIB object already exists: cannot set parameters!')
            if potential is not None:
                args[0].potential = potential
            return args[0]

        elif isinstance(args[0], str):

            # Check if the first part of the argument is an executable
            exe = shlex.split(args[0])[0]
            if shutil.which(exe) is not None:
                if lammps_command is not None:
                    raise ValueError('lammps_command seems to have been given twice')
                lammps_command = args[0]
                        
            else:
                if name != '':
                    raise ValueError('name seems to have been given twice')
                name = args[0]
        elif isinstance(args[0], lammps.lammps):
            raise TypeError('use an atomman.lammps.LAMMPSLIB object instead of a lammps.lammps object!')
        else:
            raise TypeError('unsupported argument type: should be str, LAMMPSEXE or LAMMPSLIB')

    elif len(args) > 1:
        raise ValueError('this function can only interpret one non-keyword argument')

    # Build a LAMMPSEXE object
    if lammps_command is not None:
        if name != '' or cmdargs is not None or ptr is not None or comm is not None:
            raise ValueError('name, cmdargs, ptr, and comm cannot be used with lammps_command!')
        return LAMMPSEXE(lammps_command=lammps_command, mpi_command=mpi_command, potential=potential)
    
    if mpi_command is not None:
        raise ValueError('use comm to set up MPI LAMMPS library jobs instead of mpi_command!')

    # Build a LAMMPSLIB object
    return LAMMPSLIB(name=name, cmdargs=cmdargs, ptr=ptr, comm=comm,
                     potential=potential)


