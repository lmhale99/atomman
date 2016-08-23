import atomman as am
import atomman.unitconvert as uc
import numpy as np

try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False
    
def load(pmg_struct):
    """Convert a pymatgen.Structure into an atomman.System and list of elements."""
    
    assert has_pmg, 'pymatgen not imported'

    box = am.Box(vects = pmg_struct.lattice.matrix)
    atoms = am.Atoms(natoms=pmg_struct.num_sites, prop=pmg_struct.site_properties)
    atoms.prop(key = 'pos', value = pmg_struct.cart_coords)        

    all_elements =  np.array([str(symbol) for symbol in pmg_struct.species])
    elements, atype = np.unique(all_elements, return_inverse = True)
    atype += 1
    
    atoms.prop(key='atype', value=atype)
    
    return am.System(box=box, atoms=atoms), elements
    
def dump(system, elements):
    """Convert an atomman.System and list of elements into a pymatgen.Structure."""
    
    assert has_pmg, 'pymatgen not imported'

    elements = np.asarray(elements)
    atype = system.atoms_prop(key='atype')
    
    latt = pmg.Lattice(system.box.vects)
    species = elements[atype-1]
    sites = system.atoms_prop(key='pos', scale=True)
    
    site_properties = {}
    for prop in system.atoms_prop():
        if prop != 'atype' and prop != 'pos':
            site_properties[prop] = system.atoms_prop(key=prop)

    return pmg.Structure(latt, species, sites, site_properties=site_properties)
