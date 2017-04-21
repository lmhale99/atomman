import os
import glob
import shutil
import subprocess as sp

from .Log import Log


def run(lammps_command, script_name, mpi_command=None, restart_script_name=None, return_style='object', logfile='log.lammps', flatten=None):
    """
    Calls LAMMPS to run. Returns data model containing LAMMPS output.
    
    Arguments:
    lammps_command -- The LAMMPS inline run command (sans -in script_name).
    script_name -- Path of the LAMMPS input script to use.
    
    Keyword Arguments:
    mpi_command -- The MPI inline command to run LAMMPS in parallel. Default is
                   None (no mpi).
    restart_script_name -- alternative script to use for restarting if logfile 
                           already exists. Default is None (no restarting)
    return_style -- format for the returned data. Default value is 'object'.
        'object' -- returns an atomman.lammps.Log object
        'model' -- returns a DataModelDict
    logfile -- specifies the path to the logfile to write to. Default value is
               'log.lammps'
    flatten -- indicates if the simulation data should be flattened to a 
               single table. Default value is None (no flattening)
        None -- leave all individual runs/restarts as separate tables
        'first' -- only use the first entry for a given Step value
        'last' -- only use the last entry for a given Step value
    
    lammps_command, script_name, and mpi_command can be single strings
    or lists of strings broken where spaces occur.
    """
    
    #check if restart_script_name is given
    if restart_script_name is not None:

        #check if simulation was previously started by looking for log.lammps
        if os.path.isfile(logfile):
            logname, logext = os.path.splitext(logfile)
            #Replace script_name with restart_script_name
            script_name = restart_script_name
            
            #Search for any earlier log files with the name log-*.lammps
            logids = []
            for oldlog in glob.iglob(logname + '-*' + logext):
                logids.append(int(os.path.splitext(os.path.basename(oldlog))[0][len(logname)+1:]))
            
            #Rename old logfile to keep it from being overwritten
            if len(logids) == 0:
                lognum = 1
            else:
                lognum = max(logids)+1
            shutil.move(logfile, logname + '-' + str(lognum) + logext)
        else:
            lognum = 0
    else:
        lognum = 0

    #convert lammps_command into list of terms
    if isinstance(lammps_command, (str, unicode)):
        lammps_command = lammps_command.split(' ')
    elif not isinstance(lammps_command, list):
        lammps_command = [lammps_command]   
        
    #convert script_name into list of terms
    if isinstance(script_name, (str, unicode)):
        script_name = script_name.split(' ')
    elif not isinstance(script_name, list):
        script_name = [script_name]
    
    #convert mpi_command into list of terms
    if mpi_command is None:
        mpi_command = []
    elif isinstance(mpi_command, (str, unicode)):
        mpi_command = mpi_command.split(' ')
    elif not isinstance(mpi_command, list):
        mpi_command = [mpi_command]
        
    #extra terms
    extra = []
        
    #check if logfile is not log.lammps
    if logfile != 'log.lammps':
        #convert logfile into list of terms
        extra += ['-log'] + logfile.split(' ')
        
    #try to run lammps as a subprocess
    try:
        output = sp.check_output(mpi_command + lammps_command + extra + ['-in'] + script_name)
    
    #pass LAMMPS error to Python error if failed
    except sp.CalledProcessError as e:
        if e.output != '':
            lines = e.output.split('\n')
            raise ValueError('Invalid LAMMPS input: \n%s' % lines[-2])
        else:
            raise OSError('Failed to run LAMMPS')
    
    #initialize Log object
    log = Log()
    
    #Read in from all old log files
    for i in xrange(1, lognum+1):
        log.read(logname + '-' + str(i) + logext)
        
    #Read in from current logfile
    log.read(logfile)
    
    if flatten is not None:
        log.flatten(flatten)

    if return_style == 'object':
        return log
    elif return_style == 'model':
        return log.model()
    else:
        raise ValueError('Invalid return_style')
    
    
