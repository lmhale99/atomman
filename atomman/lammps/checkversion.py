import datetime
import shlex
from pathlib import Path

# atomman imports
from . import run, run_libtest
from ..typing import lammps, lammpsapp

def checkversion(lammps_command: lammpsapp) -> dict:
    """
    Gets the LAMMPS version and date information by testing lammps_command.
    
    Parameters
    ----------
    lammps_command: str or lammps.lammps
        A LAMMPS executable or LAMMPS Python interface.
        
    Returns
    -------
    version_info: dict
        Dictionary containing 'version', the str LAMMPS version, and
        'date', the corresponding  datetime.date for the LAMMPS 
        version.
    
    Raises
    ------
    ValueError
        If lammps fails to run
    """
    if isinstance(lammps_command, lammps.lammps):
        log = run_libtest()
    else:
        # Switch for LAMMPSEXE objects
        if hasattr(lammps_command, 'lammps_command'):
            lammps_command = lammps_command.lammps_command
        
        # Use shlex to strip away any command line arguments
        lammps_command = Path(shlex.split(lammps_command)[0]).as_posix()
        
        # Run lammps_command with empty script and no log file
        tries = 0
        success = False
        while tries < 10:
            try:
                log = run(lammps_command, script='', logfile=None)
            except OSError as err:
                tries += 1
                if tries >= 10:
                    raise ValueError(f'Failed to run simulation with lammps_command {lammps_command}') from err
            except:
                raise ValueError(f'Failed to run simulation with lammps_command {lammps_command}')
            else:
                success = True
                break
        if not success:
            raise ValueError(f'Failed to run simulation with lammps_command {lammps_command} after 10 tries')
        
    # Extract lammps version and date info
    version_info = {}
    version_info['version'] = log.lammps_version
    version_info['date'] = log.lammps_date
    return version_info

def versiondate(lmp: lammpsapp) -> datetime.date:
    """
    Gets the LAMMPS version and date information from a lammps object or by
    testing a lammps_command.
    
    Parameters
    ----------
    lmp: str or lammps
        A LAMMPS executable or a lammps object.
        
    Returns
    -------
    datetime.date
        The datetime.date for the LAMMPS version.
    """
    if isinstance(lmp, lammps.lammps):
        return datetime.date.fromisoformat(str(lmp.version()))
    elif isinstance(lmp, str):
        return checkversion(lmp)['date']
    else:
        raise ValueError('lmp must be a lammps object or a lammps command string')