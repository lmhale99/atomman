# Standard Python libraries
import sys
from pathlib import Path
import shlex
import subprocess
from typing import Optional, Tuple

# atomman imports
from . import Log, LammpsError

def run(lammps_command: str,
        script_name: Optional[str] = None,
        script: Optional[str] = None,
        mpi_command: Optional[str] = None,
        restart_script_name: Optional[str] = None,
        restart_script: Optional[str] = None,
        restart_filename: Optional[str] = None,
        restart_lognum: Optional[int] = None,
        partition: Optional[str] = None,
        logfile: Optional[str] = 'log.lammps',
        screen: bool = True,
        return_log: bool = True,
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
        Restart scripts will be used instead of the regular scripts if they are
        given, logfile already exists, and if restart_filename is given at least
        one matching file is found.  Requires logfile to not be None.
    restart_script : str, optional
        Alternate LAMMPS script command lines to use for restart runs.
        Restart scripts will be used instead of the regular scripts if they are
        given, logfile already exists, and if restart_filename is given at least
        one matching file is found.  Requires logfile to not be None.
    restart_filename : str or None, optional
        File name (with wildcards) indicating the restart files the simulation
        would start from if they exist.  If given, restart scripts will not be
        used unless at least one matching file is found to exist.
    restart_lognum : int or None, optional
        Handling of restart scripts is done in two steps: initial checks
        for existing log files and reading in all log files after running.
        If the initial checks are handled externally beforehand, the lognum
        setting (i.e. how many previous log files there are) can be passed in
        allowing for them to still all be read in afterwards.  Note: this cannot
        be used with restart_script or restart_script_name so it is your
        responsibility to set script or script_name based on if it is a
        restarted calculation or not!
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
    return_log : bool, optional
        Setting this to False skips reading in any log files and makes the
        method return None.  Useful for complex simulations in which you
        do not analyze the raw log data.
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
        if restart_lognum is not None:
            raise  ValueError('cannot give restart_lognum with restart_script or restart_script_name: restart is assumed!')
        if restart_script_name is not None and restart_script is not None:
            raise  ValueError('Cannot give both restart_script and restart_script_name')
        if logfile is None:
            raise ValueError('logfile must be given to automatically determine restart status')
        
        # Count previous log files
        logname, logext, lognum = restart_check(logfile, restart_filename)

        if lognum > 0:
            # Replace script parameters with restart_script parameters
            script = restart_script
            script_name = restart_script_name

    # Set lognum value if given
    elif restart_lognum is not None:
        lognum = restart_lognum
        logname = logext = None

    # Set default lognum value
    else:
        lognum = 0
        logname = logext = None
    
    # Initialize the run command
    command = ''

    # Add mpi_command
    if mpi_command is not None:
        command += mpi_command + ' '
    
    # Add lammps_command
    terms = shlex.split(lammps_command)
    terms[0] = Path(terms[0]).as_posix()
    command += ' '.join(terms) + ' '

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
        if screen is False and logfile != 'none':
            with open(logfile) as f:
                output = f.read()
        else:
            output = e.output
        raise LammpsError(output) from e
    
    # Read log contents from screen output if generated
    if screen is True:
        logfile = output.stdout

    if return_log:
        # Read in all log files to a Log object
        log = read_logs(logfile, logname, logext, lognum)
    
        return log

def restart_check(logfile,
                  restart_filename: Optional[str] = None
                  ) -> Tuple[str, str, int]:
    """
    Checks for existing log and restart files to determine if a previously
    incomplete simulation is to be restarted.  If found, the previous log file
    will be renamed to prevent it from being overwritten.

    USAGE: This is used internally by run().  If you wish to use it separately,
    pass the returned lognum value to run() or read_logs().

    Parameters
    ----------
    logfile : str
        The path to the logfile to check for.
    restart_filename : str or None
        If given, the check will also verify that at least one matching restart
        files also exists.  Ideally, this should always be given but is optional
        for backwards compatibility with older atomman versions.
    
    Returns
    -------
    logname : str
        The log file prefix of the previous simulation runs, which matches the
        given logfile's stem.
    logext : str
        The log file extension of the previous simulation runs, which matches the
        given logfile's extension.
    lognum : int
        The number of previous log files found.  If 0, then a new simulation
        should be started rather than a restart.
    """

    # Count existing restart files
    if restart_filename is None:
        num_restarts = 1
    else:
        num_restarts = len([fname for fname in Path().glob(restart_filename)])

    # Check if simulation was previously started by looking for the logfile
    # (and at least one restart file)
    logfile = Path(logfile)
    if logfile.is_file() and num_restarts > 0:
        logname = logfile.stem
        logext = logfile.suffix
        
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
        logname = logext = None

    return logname, logext, lognum

def read_logs(logfile: str,
              logname: Optional[str] = None,
              logext: Optional[str] = None,
              lognum: int = 0) -> Log:
    """
    Creates a Log object based on either a single LAMMPS run or multiple
    LAMMPS restart runs.  For single runs, this is essentially equivalent
    to Log.read().

    Parameters
    ----------
    logfile : str
        The log file name or contents for the final LAMMPS run.
    logname : str or None, optional
        The log file prefix that any previous log files have.
        Previous log files should be named {logname}-{i}{logext},
        where i ranges from 1 up to and including lognum.  Optional
        if logfile is a file named {logname}{logext}.
    logext : str or None, optional
        The log file extension that any previous log files have.
        Previous log files should be named {logname}-{i}{logext},
        where i ranges from 1 up to and including lognum.  Optional
        if logfile is a file named {logname}{logext}.
    lognum : int, optional
        The total number of previous log files.
        Previous log files should be named {logname}-{i}{logext},
        where i ranges from 1 up to and including lognum.  Default
        value is 0, i.e. no previous log files.
    
    Returns
    -------
    atomman.lammps.Log
        The interpreted log contents for the log file(s).
    """
    # Initialize Log object
    log = Log()

    # Identify logname and logext if needed
    if lognum > 0 and (logname is None or logext is None):
        logfile = Path(logfile)
        if not logfile.is_file():
            raise ValueError('logname and logext must be given if logfile is not a file name!')
        logname = logfile.stem
        logext = logfile.suffix
    
    # Read in all old runs
    for i in range(1, lognum+1):
        log.read(f'{logname}-{i}{logext}')
    
    # Read in current run
    log.read(logfile)
    
    return log


def run_libtest() -> Log:
    """
    Runs a subprocess to test building a new lammps.lammps object.  Used to
    verify the LAMMPS Python interface works and extract version information
    """

    # Build Python script that creates a lammps.lammps object and does no simulation.
    script = '\n'.join(["import lammps", "lmp = lammps.lammps(cmdargs=['-l', 'none'])", "lmp.close()"])
    output = subprocess.run([sys.executable], input=script, check=True, capture_output=True, text=True)
    
    # Initialize Log object and read screen output
    log = Log()
    log.read(output.stdout)
    
    return log
    