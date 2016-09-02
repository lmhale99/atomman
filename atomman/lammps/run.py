import subprocess as sp
from log_extract import log_extract

def run(lammps_command, script_name, mpi_command=None):
    """
    Calls LAMMPS to run. Returns data model containing LAMMPS output.
    
    Arguments:
    lammps_command -- The LAMMPS inline run command (sans -in script_name).
    script_name -- Path of the LAMMPS input script to use.
    mpi_command -- (optional) The MPI inline command to run LAMMPS in parallel.
    
    All of the arguments can be either strings or lists of strings. 
    If a term contains spaces, then either leave the spaces in (no 
    additional outer quotes) or split into a list based on the spaces.
    """
    
    #convert lammps_command into list of terms
    if isinstance(lammps_command, (str, unicode)) and ' ' in lammps_command:
        lammps_command = lammps_command.split(' ')
    elif not isinstance(lammps_command, list):
        lammps_command = [lammps_command]   
        
    #convert script_name into list of terms
    if isinstance(script_name, (str, unicode)) and ' ' in script_name:
        script_name = script_name.split(' ')
    elif not isinstance(script_name, list):
        script_name = [script_name]
    
    #convert mpi_command into list of terms
    if mpi_command is None:
        mpi_command = []
    elif isinstance(mpi_command, (str, unicode)) and ' ' in mpi_command:
        mpi_command = mpi_command.split(' ')
    elif not isinstance(mpi_command, list):
        mpi_command = [mpi_command]
    
    #try to run lammps as a subprocess and return log_extract output
    try:
        return log_extract(sp.check_output(mpi_command + lammps_command + ['-in'] + script_name))
    
    #pass LAMMPS error to Python error if failed
    except sp.CalledProcessError as e:        
        if e.output != '':
            lines = e.output.split('\n')
            raise ValueError('Invalid LAMMPS input: \n%s' % lines[-2])
        else:
            raise OSError('Failed to run LAMMPS')
