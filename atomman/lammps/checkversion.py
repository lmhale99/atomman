# coding: utf-8

# atomman imports
from . import run

def checkversion(lammps_command):
    """
    Gets the LAMMPS version and date information by testing lammps_command.
    
    Parameters
    ----------
    lammps_command: str
        A LAMMPS executable.
        
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
    # Run lammps_command with empty script and no log file
    try:
        log = run(lammps_command, script='', logfile=None)
    except:
        raise ValueError(f'Failed to run simulation with lammps_command {lammps_command}')
        
    # Extract lammps version and date info
    version_info = {}
    version_info['version'] = log.lammps_version
    version_info['date'] = log.lammps_date
    return version_info