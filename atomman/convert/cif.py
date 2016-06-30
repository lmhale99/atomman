import atomman as am
import numpy as np

try:
    import diffpy.Structure
    has_diffpy = True
except:
    has_diffpy = False

def load(cif):
    """Reads in a CIF crystal file and returns an atomman.System and list of elements."""
    
    assert has_diffpy, 'diffpy.Structure not imported'
    
    dps = diffpy.Structure.structure.Structure()
    
    if hasattr(cif, 'read'):
        cif = cif.read()
        dps.readStr(cif)
    else:
        dps.read(cif)
    
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