import numpy as np

import atomman.unitconvert as uc
from ...lammps.style import unit as get_lammps_units

def dump(
    system,
    units: str = 'metal',
    include_velocities: bool = True) -> dict:
    """
    Extracts data stored in an atomman System and transforms it into simple
    terms that can be directly passed to the pair_lib_info() methods of
    PotentialLAMMPS and PotentialLAMMPSKIM objects.

    Typical users will want to use the lammps_lib dump style instead as it
    is simpler.  This is used by iprPy_fit for a minor speedup by extracting
    system parameters before evaluating the same systems repeatedly.
    
    Parameters
    ----------
    system : atomman.System
        The system that LAMMPS will replicate.
    units : str, optional
        The LAMMPS units option associated with the data file.  If neither
        units or potential is given, will set units 'metal'.
    include_velocities : bool, optional
        Indicates if velocity information in the system (if present) is to
        be extracted and included in the returned outputs.  Default value is
        True.
    
    Returns
    -------
    params : dict
        Dict containing the following by name.
    pbc : list
        Tuple of three bool indicating which directions are periodic.
    region_params : list
        The LAMMPS parameters to define the box shape and position:
        [xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz].
    atype : list
        The atype of each atom.
    x : list
        The flattened array of atom coordinates.
    v : list or None
        The flattened array of atom velocities. Will be None if
        include_velocities is False.
    """
    lammps_units = get_lammps_units(units)

    # Test that box parameters are compatible with LAMMPS
    if not system.box.is_lammps_norm():
        raise ValueError('System not normalized for LAMMPS compatibility')
    
    region_params = []
    for param in ['xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi', 'xy', 'xz', 'yz']:
        region_params.append(uc.get_in_units(getattr(system.box, param), lammps_units['length']))
    
    # Convert pbc to LAMMPS format
    bflags = np.array(['m','m','m'])
    bflags[system.pbc] = 'p'
    pbc = bflags
    
    atype = system.atoms.atype
    x = uc.get_in_units(system.atoms.pos, lammps_units['length']).flatten()

    if include_velocities and 'velocity' in system.atoms.prop():
        v = uc.get_in_units(system.atoms.velocity, lammps_units['velocity']).flatten()
    else:
        v = None

    return dict(
        pbc = pbc,
        natoms = system.natoms,
        natypes = system.natypes,
        region_params = region_params,
        atype = atype,
        symbols=system.symbols,
        x = x,
        v = v
    )
