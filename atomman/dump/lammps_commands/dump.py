# coding: utf-8

# Standard Python libraries
from copy import deepcopy
from typing import Optional

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from potentials.record.BasePotentialLAMMPS import BasePotentialLAMMPS

def dump(system,
         atom_style: Optional[str] = None,
         units: Optional[str] = None,
         natypes: Optional[int] = None, 
         potential: Optional[BasePotentialLAMMPS] = None,
         uvws: Optional[npt.ArrayLike] = None,
         shift: Optional[npt.ArrayLike] = None,
         size_mults: Optional[npt.ArrayLike] = None,
         return_pair_info: bool = False,
         safecopy: bool = False) -> str:
    """
    Write LAMMPS commands that result in LAMMPS generating an atomic configuration
    based on the supplied system.  
    
    NOTE: It is recommended that this only be used for small systems (e.g. unit cells)
    and to instead use  'atom_data' dump style for more complex configurations.
    
    Parameters
    ----------
    system : atomman.System 
        The system that LAMMPS will replicate.
    atom_style : str, optional
        The LAMMPS atom_style option associated with the data file.  If neither
        atom_style or potential is given, will set atom_style to 'atomic'.
    units : str, optional
        The LAMMPS units option associated with the data file.  If neither
        units or potential is given, will set units 'metal'.
    natypes : int, optional
        Allows the natypes value to be manually changed.  This is needed if
        natypes needs to be greater than the current number of atypes.  If
        neither natypes or potential is given, will use system.natypes.
    potential : atomman.lammps.Potential, optional
        Potential-specific values of atom_style, units, and natypes can be
        extracted from a Potential object.  If both potential and any of the
        individual values are given, the individual values will be used.
    uvws : array-like object, optional
        3x3 array of integer Miller crystal vectors to reorient the system along.
        Currently only supported for cubic systems.
    shift : array-like object, optional
        A 3D rigid-body shift vector to apply to the atomic positions.  Useful to
        avoid positioning atoms directly on the boundaries.  Default value is 
        [0.01, 0.01, 0.01].
    size_mults : array-like object, optional
        3 sets of (-,+) size multipliers to apply to the system.  If not given, no
        multiplications of the system will be made.
    return_pair_info : bool, optional
        Indicates if the LAMMPS command lines associated with setting mass,
        pair_style and pair_coeff are included in the returned info.  If True,
        potential must be given.  Default value is False.
    safecopy : bool, optional
        The LAMMPS commands requires all atoms to be inside box bounds, i.e.
        "wrapped".  If safecopy is True then a copy of the system is made to
        keep the original unwrapped.  Default value is False.
    
    Returns
    -------
    str
        The LAMMPS input command lines to generate the system for itself.
        If return_pair_info is True and potential is given, the LAMMPS input
        command lines for the potential will also be included.

    Raises
    ------
    ValueError
        If return_pair_info is True and potential is not given.
    ValueError
        If the given system is not normalized for LAMMPS compatibility.
    ValueError
        If rotation axes are given for a non-orthogonal system.
    """
    
    # Test that box parameters are compatible with LAMMPS
    if not system.box.is_lammps_norm():
        raise ValueError('System not normalized for LAMMPS compatibility')
    
    # Add comment line
    info = '# Script prepared using atomman Python package\n\n'
    
    # Extract potential-based parameters
    if potential is not None:
        if units is None:
            units = potential.units
        if atom_style is None:
            atom_style = potential.atom_style
        if natypes is None:
            natypes = len(potential.normalize_symbols(system.symbols))
    
    # Set default parameter values
    else:
        if units is None:
            units = 'metal'
        if atom_style is None:
            atom_style = 'atomic'
        if natypes is None:
            natypes = system.natypes
    
    # Generate units and atom_style lines
    info += f'units {units}\n'
    info += f'atom_style {atom_style}\n\n'
    
    # Wrap atoms
    if safecopy:
        system = deepcopy(system)
    system.wrap()
    
    family = system.box.identifyfamily()
    
    # Generate boundary condition
    info += 'boundary '
    for pbc in system.pbc:
        if pbc:
            info += 'p '
        else:
            info += 'm '
    info += '\n\n'
    
    # Start lattice command with box parameters
    vects = system.box.vects
    info += 'lattice custom 1.0 &\n'
    info += f'        a1 {vects[0,0]:.12f} {vects[0,1]:.12f} {vects[0,2]:.12f} &\n'
    info += f'        a2 {vects[1,0]:.12f} {vects[1,1]:.12f} {vects[1,2]:.12f} &\n'
    info += f'        a3 {vects[2,0]:.12f} {vects[2,1]:.12f} {vects[2,2]:.12f} &\n'
    
    # Add origin shift
    if shift is None:
        shift = np.array([0.01, 0.01, 0.01])
    info += f'        origin {shift[0]} {shift[1]} {shift[2]} &\n'
    
    # Add spacing and orient
    if uvws is not None:
        if family not in ['cubic', 'trigonal', 'orthorhombic']:
            raise ValueError('uvws not currently supported for non-orthorhombic systems')
         
        # Adjust crystal spacing for systems to be (nearly) perfectly periodic across boundaries
        spacing = np.zeros(3)
        for i in range(3):
            spacing[i] = vects[i,i] * np.linalg.norm(uvws[i])
            
        info += f'        spacing {spacing[0]:.12f} {spacing[1]:.12f} {spacing[2]:.12f}&\n'
        info += f'        orient x {uvws[0,0]} {uvws[0,1]} {uvws[0,2]} &\n'
        info += f'        orient y {uvws[1,0]} {uvws[1,1]} {uvws[1,2]} &\n'
        info += f'        orient z {uvws[2,0]} {uvws[2,1]} {uvws[2,2]} &\n'
    
    # Add basis (i.e. positions)
    spos = system.atoms_prop('pos', scale=True)
    for i in range(system.natoms):
        info += f'        basis {spos[i,0]:.12f} {spos[i,1]:.12f} {spos[i,2]:.12f}'
        if i < system.natoms - 1:
            info += ' &\n'
    info += '\n\n'
    
    # Generate region
    # Identify if system is orthorhombic or not
    if family in ['cubic', 'trigonal', 'orthorhombic']:
        region_style = 'block'
    else:
        region_style = 'prism'
    
    # Generate the system region 
    if size_mults is None:
        size_mults = np.array([[0,1], [0,1], [0,1]])
    else:
        size_mults = np.asarray(size_mults)
    info += f'region box {region_style} {size_mults[0,0]} {size_mults[0,1]} {size_mults[1,0]} {size_mults[1,1]} {size_mults[2,0]} {size_mults[2,1]}'
    if region_style == 'prism':
        size_xy = vects[1][0] * (size_mults[0,1] - size_mults[0,0]) / system.box.a
        size_xz = vects[2][0] * (size_mults[0,1] - size_mults[0,0]) / system.box.a
        size_yz = vects[2][1] * (size_mults[1,1] - size_mults[1,0]) / system.box.b
        info += f' {size_xy} {size_xz} {size_yz}'
    info += '\n'
    
    # Create box and atoms
    info += f'create_box {natypes} box\n'
    info += 'create_atoms 1 box '
    for i in range(system.natoms):
        if system.atoms.atype[i] > 1:
            info += f' &\n             basis {i+1} {system.atoms.atype[i]}'
    info += '\n'
    
    # Set pair_info
    if return_pair_info is True:
        if potential is None:
            raise ValueError('return_pair_info = True requires that potential be given')
        
        info += '\n'
        info += potential.pair_info(symbols=system.symbols,
                                    masses=system.masses)
    
    return info