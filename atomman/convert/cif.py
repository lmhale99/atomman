import diffpy.Structure 
import atomman as am
import numpy as np


def load(cif):
    if hasattr(cif, 'read'):
        cif = cif.read()
    
    dps = diffpy.Structure.structure.Structure()
    dps.readStr(cif)
    
    all_elements = dps.element
    elements, all_atype = np.unique(all_elements, return_inverse=True)
    all_atype += 1
    all_pos = dps.xyz
    
    atype = []
    pos = []
    for a1, p1 in zip(all_atype, all_pos):
        noMatch = True
        for a2, p2 in zip(atype, pos):
            if a1 == a2 and np.allclose(p1, p2):
                noMatch = False
                break
        if noMatch:
            atype.append(a1)
            pos.append(p1)
    
    atoms = am.Atoms(natoms=len(pos), prop={'atype':atype, 'pos':pos})    
    
    box = am.Box(a=dps.lattice.a, 
                 b=dps.lattice.b, 
                 c=dps.lattice.c, 
                 alpha=dps.lattice.alpha, 
                 beta=dps.lattice.beta, 
                 gamma=dps.lattice.gamma)
    
    return am.System(atoms=atoms, box=box, pbc=(True,True,True), scale=True), list(elements)
   
    
if __name__ == '__main__':
    import glob
    import os
    for fname in glob.iglob('E:\\crystal-structures\\oqmd\\cif\\*.cif'):
        with open(fname) as f:
            elements, system = load(f)
    
        print os.path.basename(fname), system.box.a, system.box.b, system.box.c, system.box.alpha, system.box.beta, system.box.gamma