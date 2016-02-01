from subprocess import check_output
from log_extract import log_extract

def run(lammps_exe, script_name, mpi_command=None):
    """
    Calls LAMMPS to run. Returns data model containing LAMMPS output.
    
    Arguments:
    lammps_exe -- name (and location) of the LAMMPS executable to use.
    script_name -- name of the LAMMPS input script to use.
    mpi_command -- (optional) command line associated with calling MPI to run LAMMPS.
    """
    if mpi_command is None:
        return log_extract(check_output(lammps_exe + ' -in ' + script_name, shell=True))
    else:
        return log_extract(check_output(mpi_command + ' ' + lammps_exe + ' -in ' + script_name, shell=True))
    