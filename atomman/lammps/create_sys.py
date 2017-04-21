import os
import numpy as np
from run import run
import atom_dump
   
def create_sys(lammps_command, system_info):
    """
    Uses LAMMPS to generate a System based on the supplied system_info.
    
    Arguments:
    lammps_command -- the lammps command/executable to use.
    system_info -- LAMMPS input script command lines associated with creating a new system.
    
    system_info can be generated using atomman.lammps.sys_gen.
    """

    newline = '\n'
    script = newline.join([system_info,
                           '',
                           'mass * 1',
                           'pair_style none',
                           'atom_modify sort 0 0.0',
                           '',                           
                           'dump dumpit all custom 100000 temp.dump id type x y z',
                           'run 0'])
    f = open('create_sys.in', 'w')
    f.write(script)
    f.close()
    
    output = run(lammps_command, 'create_sys.in')

    system = atom_dump.load('temp.dump')
    os.remove('create_sys.in')
    os.remove('log.lammps')
    os.remove('temp.dump')
    os.remove('temp.dump.json')
    return system      