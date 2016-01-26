from subprocess import check_output
from log_extract import log_extract

def run(lammps_exe, script_name, mpi_command=None):
    if mpi_command is None:
        return log_extract(check_output(lammps_exe + ' -in ' + script_name, shell=True))
    else:
        return log_extract(check_output(mpi_command + ' ' + lammps_exe + ' -in ' + script_name, shell=True))
    