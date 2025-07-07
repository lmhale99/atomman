# coding: utf-8
# Standard Python libraries
from pathlib import Path
import shlex
import subprocess
from typing import Optional

# atomman imports
from . import Log, LammpsError

def run(lammps_command: str,
        script_name: Optional[str] = None,
        script: Optional[str] = None,
        mpi_command: Optional[str] = None,
        restart_script_name: Optional[str] = None,
        restart_script: Optional[str] = None,
        partition: Optional[str] = None,
        logfile: str = 'log.lammps',
        screen: bool = True,
        suffix: Optional[str] = None) -> Log:
    """
    Calls LAMMPS to run. Returns a Log object of the screen/logfile output.
    
    Parameters
    ----------
    lammps_command : str
        The LAMMPS inline run command (sans -in script_name).
    script_name : str, optional
        Path of the LAMMPS input script file to use.  Either script_name or
        script must be given.
    script : str, optional
        The LAMMPS input script command lines to use.  Either script_name or
        script must be given.
    mpi_command : str, optional
        The MPI inline command to run LAMMPS in parallel. Default value is 
        None (no mpi).
    restart_script_name : str, optional
        Path to an alternate LAMMPS input script file to use for restart runs.
        If given, the restart script will be used if the specified logfile
        already exists.  Requires logfile to not be None.
    restart_script : str, optional
        Alternate LAMMPS script command lines to use for restart runs.
        If given, the restart script will be used if the specified logfile
        already exists.  Requires logfile to not be None.
    partition : str or None, optional
        The LAMMPS partition setting to use for the calculation.  This is
        required for calculations like NEB that run multiple simulations at
        the same time.
    logfile : str or None, optional
        Specifies the path to the logfile to write to.  Default value is
        'log.lammps'.  If set to None, then no logfile will be created.
    screen : bool, optional
        If True (default), then the resulting Log object is built from the
        LAMMPS screen output.  If False, then LAMMPS outputs no screen info
        and the Log object will be built by reading logfile.
    suffix : str, optional
        Allows for the LAMMPS suffix option to be specified to use any of the
        accelerated versions of pair styles if available.
    
    Returns
    -------
    atomman.lammps.Log 
        Contains the processed screen/logfile contents.  Will not be returned
        if logfile is None and screen is False as there would be no means to
        capture the LAMMPS output.
    """

    # Check if either restart_script_name or restart_script is given
    if restart_script_name is not None or restart_script is not None:
        if restart_script_name is not None and restart_script is not None:
            raise  ValueError('Cannot give both restart_script and restart_script_name')
        if logfile is None:
            raise ValueError('logfile must be given to automatically determine restart status')
        else:
            logfile = Path(logfile)
        
        # Check if simulation was previously started by looking for the logfile
        if logfile.is_file():
            logname = logfile.stem
            logext = logfile.suffix
            
            # Replace script parameters with restart_script parameters
            script = restart_script
            script_name = restart_script_name
            
            # Search for any earlier log files with the name log-*.lammps
            maxlogid = 0
            for oldlog in Path().glob(f'{logname}-*{logext}'):
                logid = int(oldlog.stem.split('-')[-1])
                if logid > maxlogid:
                    maxlogid = logid
            
            # Rename old logfile to keep it from being overwritten
            lognum = maxlogid + 1
            logfile.rename(f'{logname}-{lognum}{logext}')
        else:
            lognum = 0
    else:
        lognum = 0
    

    # Initialize the run command
    command = ''

    # Add mpi_command
    if mpi_command is not None:
        command += mpi_command + ' '
    
    # Add lammps_command
    command += f'"{Path(lammps_command).as_posix()}" '

    # Add logfile
    if logfile is None:
        logfile = 'none'
    if logfile != 'log.lammps':
        command += f'-log {logfile} '

    # Add script_name
    if script_name is not None:
        if script is not None:
            raise  ValueError('Cannot give both script and script_name')
        command += f'-in {script_name} '
    elif script is None:
        raise ValueError('script or script_name must be given')

    # Add partition
    if partition is not None:
        command += f'-partition {partition} '

    # Add screen
    if screen is False:
        command += f'-screen none '

    # Add suffix
    if suffix is not None:
        command += f'-suffix {suffix} '

    # Use shlex to split command
    command = shlex.split(command)
    
    # Try to run lammps as a subprocess
    try:
        output = subprocess.run(command, input=script, check=True,
                                capture_output=True, text=True)
    
    # Convert LAMMPS error to a Python error if failed
    except subprocess.CalledProcessError as e:
        raise LammpsError(e.output) from e
    
    # Initialize Log object
    log = Log()
    
    # Read in all old runs
    for i in range(1, lognum+1):
        log.read(f'{logname}-{i}{logext}')
    
    # Read in current run
    if screen:
        log.read(output.stdout)
    else:
        log.read(logfile)
    
    return log