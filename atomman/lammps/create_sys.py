import os
import numpy as np
from atomman.lammps import atom_dump, run
   
def create_sys(lammps_exe, system_info):
    """
    Uses LAMMPS to generate a System based on the supplied system_info.
    
    Arguments:
    lammps_exe -- name (and location) of the LAMMPS executable to use.
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
                           'dump_modify dumpit format "%i %i %.13e %.13e %.13e"',                           
                           'run 0'])
    f = open('create_sys.in', 'w')
    f.write(script)
    f.close()
    output = run(lammps_exe, 'create_sys.in')

    system = atom_dump.load('temp.dump')
    os.remove('create_sys.in')
    os.remove('log.lammps')
    os.remove('temp.dump')
    return system      