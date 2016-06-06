from copy import deepcopy
import numpy as np

def normalize(system):
    """
    The normalize function takes any arbitrary system and transforms it to
    be comptatible with LAMMPS.  In particular, LAMMPS systems must have:
        1. Right-handed box vectors.
        2. avect = [lx, 0.0, 0.0]
        3. bvect = [xy, ly,  0.0]
        4. cvect = [xz, yz,  lz]
        5. All atoms initially inside the box dimensions.
    Note: large box tilt factors are not adjusted with this function.
    As such, the LAMMPS command 'box tilt large' may be needed.
    """
    
    #create a copy of the system 
    system = deepcopy(system)
        
    #Swap z-axis direction if box is left-handed.
    if np.dot(np.cross(system.box.avect, system.box.bvect), system.box.cvect) < 0:
        spos = system.atoms_prop(key='pos', scale=True)
        system.box_set(avect=system.box.avect, 
                       bvect=system.box.bvect, 
                       cvect=-system.box.cvect)
        spos = np.dot(spos, [[1,0,0], [0,1,0], [0,0,-1]])
        system.atoms_prop(key='pos', value=spos, scale=True)

    #Rebuild box using a, b, c, alpha, beta, gamma    
    system.box_set(a=system.box.a, b=system.box.b, c=system.box.c,
                   alpha=system.box.alpha, beta=system.box.beta, gamma=system.box.gamma,
                   scale=True)
        
    #wrap to ensure all atoms are within the box    
    system.wrap()
    
    return system