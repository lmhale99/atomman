import ase

def r0(atoms: ase.Atoms):
    """
    Computes the shortest interatomic distance for an ase.Atoms object using
    native ase actions to mimic atomman.System.r0().

    Parameters
    ----------
    atoms : ase.Atoms
        The atomic configuration to analyze. Ideally, this atomic configuration
        should be small, like a unit cell!
    """
    # Find shortest periodic lattice parameter
    try:
        box_r0 = atoms.cell.lengths()[atoms.pbc].min()
    except:
        box_r0 = None
    
    # Find shortest interatomic vector
    if len(atoms) > 1:
        
        # Use r0 for atom 0 as cutoff
        cutoff = atoms.get_distances(0, range(1, len(atoms)), mic=True).min() * 1.01
        
        # Find smallest r0 across all atoms
        atom_r0 = ase.neighborlist.neighbor_list('d', atoms, cutoff).min()
        
    else:
        atom_r0 = None

    if box_r0 is not None and atom_r0 is not None:
        if atom_r0 < box_r0:
            return atom_r0
        else:
            return box_r0
    elif box_r0 is not None:
        return box_r0
    elif atom_r0 is not None:
        return atom_r0
    else:
        raise ValueError('No atoms to compare or periodic boundaries found!')