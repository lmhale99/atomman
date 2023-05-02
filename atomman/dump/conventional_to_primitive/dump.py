# coding: utf-8
# Standard Python libraries
from typing import Optional, Union, Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from ... import Box, System
from ...tools import miller

def dump(system: System,
         setting: str = 'p',
         smallshift: Optional[npt.ArrayLike] = None,
         rtol: float = 1e-05,
         atol: float = 1e-08,
         check_basis: bool = True,
         check_family: bool = True,
         return_transform: bool = False,
         ) -> Union[System, Tuple[System, np.ndarray]]:
    """
    Transforms a conventional unit cell system of a specified Bravais space
    lattice setting into a primitive unit cell. The
    primitive_to_conventional and conventional_to_primitive dump styles are
    meant to be inverse operations, within floating point precision, to provide
    compatible primitive and conventional unit cells.

    NOTE: This dump style expects that the original starting system is a
    conventional unit cell, and only limited checks are performed to assert
    this!  Use the 'primitive_cell' dump style for a more comprehensive
    primitive unit cell identifier.

    Parameters
    ----------
    system : atomman.System
        A conventional unit cell system to find the corresponding primitive
        unit cell for.
    setting : str, optional
        The lattice setting value.  Allowed values are 'p' for primitive,
        'i' for body-centered, 'f' for face-centered, 'a', 'b', or 'c' for
        side-centered, and 't', 't1', or 't2' for hexagonal representations of
        trigonal systems.  Giving 't' will check if system is consistent with
        either 't1' or 't2'.
    smallshift : array-like object or None, optional
        A small rigid body shift to apply to the atomic positions when searching
        for which atoms are within the primitive cell.  This helps avoid
        identification issues when atoms are directly on the box boundaries.
        The default value of None will use a smallshift of [0.001, 0.001, 0.001].
    rtol : float, optional
        Relative tolerance to use for numpy.isclose.  This is used here to check
        that the conventional cell has atoms in the expected lattice positions
        for the given setting.
    atol : float, optional
        Absolute tolerance to use for numpy.isclose.  This is used here to check
        that the conventional cell has atoms in the expected lattice positions
        for the given setting.
    check_basis : bool, optional
        If True (default), a quick check will be performed on the system to see
        if it appears consistent with a Bravais space lattice with the given
        setting.  Turning this check off may be necessary for more complex
        cases, such as non-conventional cell representations and complex unit
        cells where no atoms are at the lattice site [0, 0, 0].
    check_family : bool, optional
        If True (default), then the Bravais space lattice check will include
        a check that the crystal family is consistent with a Bravais lattice
        of the given setting. For example, Bravais lattices with setting 'f'
        only exist for cubic and orthogonal cells.  This check is not done if
        either check_family or check_basis is False.  Turning this off allows
        for transformations of non-conventional cells.
    return_transform : bool, optional
        Indicates if the Cartesian transformation matrix associated with
        rotating from the conventional cell to primitive cell orientations is
        returned.  Default value is False.

    Returns
    -------
    p_ucell : atomman.System
        The primitive unit cell obtained by transforming the given conventional
        unit cell system.
    transform : numpy.ndarray
        The Cartesian transformation matrix associated with converting from the
        primitive cell orientation to the conventional cell orientation.  Only
        returned if return_transform is True.

    Raises
    ------
    ValueError
        If smallshift is not a 3D vector.
    AssertionError
        If the algorithm fails to find the expected number of atoms in the
        primitive cell.
    """

    # Handle smallshift parameter values
    if smallshift is None:
        smallshift = np.array([0.001, 0.001, 0.001])
    else:
        smallshift = np.asarray(smallshift)
        if smallshift.shape != (3, ):
            raise ValueError('smallshift must be a 3D vector')

    # Check that system is of the proper setting
    if check_basis and setting != 't':
        is_basis = check_setting_basis(system, setting=setting, rtol=rtol, atol=atol,
                                        check_family=check_family)
        if not is_basis:
            raise ValueError('system atoms do not seem to match indicated setting')
    
    # Check for generic 't' and try to identify if t1 or t2
    elif check_basis and setting == 't':
        is_t1 = check_setting_basis(system, setting='t1', rtol=rtol, atol=atol,
                                    check_family=check_family)
        is_t2 = check_setting_basis(system, setting='t2', rtol=rtol, atol=atol,
                                    check_family=check_family)
        if is_t1:
            setting = 't1'
        elif is_t2:
            setting = 't2'
        else:
            raise ValueError('system atoms do not seem to match either t1 or t2 setting')

    # 3x3x3 primitive cell for trigonal-hexagonal transformations
    # 2x2x2 primitive cell for other transformations
    if setting in ['t1', 't2', 't']:
        multip = 3
    else:
        multip = 2

    cps_uvws = miller.vector_primitive_to_conventional(multip * np.identity(3), setting=setting)
    # Construct a 2x2x2 (or 3x3x3) primitive supercell by rotating the conventional cell
    p_scell, transform = system.rotate(cps_uvws, return_transform=True)
    # Create a box for the primitive unit cell using half of all three box vects of p8_cell
    box = Box(vects = p_scell.box.vects / multip)

    # Apply the smallshift to p8_cell atomic positions to avoid boundary issues
    p_scell.atoms.pos += smallshift
    p_scell.wrap()

    # Identify atoms inside p_box
    keepindex = box.inside(p_scell.atoms.pos)

    num_expected = int(p_scell.natoms / multip ** 3)
    assert np.sum(keepindex) == num_expected, f'{np.sum(keepindex)} atoms found, {num_expected} expected'
        
    # Reverse smallshift
    p_scell.atoms.pos -= smallshift
    p_scell.wrap()

    # Create Atoms by slicing from p8_cell
    atoms = p_scell.atoms[keepindex]
    p_ucell = System(box=box, atoms=atoms, symbols=system.symbols)
    p_ucell.wrap()

    # Search for atom near periodic (0,0,0)
    dmag = np.array(p_ucell.dmag(range(p_ucell.natoms), np.zeros((p_ucell.natoms, 3))))
    if dmag.shape == ():
        dmag = np.array([dmag])
    zeroindex = np.isclose(dmag, 0.0)

    # Adjust atoms so identified atom is exactly at (0,0,0)
    if np.sum(zeroindex) == 1:
        zerocoord = p_ucell.atoms.pos[zeroindex]
        p_ucell.atoms.pos -= zerocoord

        p_ucell.wrap()

    if return_transform:
        return p_ucell, transform
    else:
        return p_ucell

def check_setting_basis(ucell: System,
                        setting: str = 'p',
                        rtol: float = 1e-05,
                        atol: float = 1e-08,
                        check_family: bool = True) -> bool:
    """
    Checks if a unit cell system is consistent with one of the standard
    conventional cell settings.  For the indicated cell setting, a search is
    performed to verify that atoms of the same type are positioned at the
    lattice site(s), and that the lattice parameters are consistent with a
    crystal family that has that setting.

    NOTE: This is only meant as a quick check and does not perform a
    comprehensive symmetry analysis of the atomic coordinates.

    Parameters
    ----------
    ucell : atomman.System
        The unit cell system to check if it appears consistent with the given
        crystal space group lattice setting.
    setting : str
        The lattice setting value.  Allowed values are 'p' for primitive,
        'i' for body-centered, 'f' for face-centered, 'a', 'b', or 'c' for
        side-centered, and 't1' or 't2' for hexagonal representations of
        trigonal systems.
    rtol : float, optional
        Relative tolerance to use for numpy.isclose.  This is used both to
        check for atoms in the equivalent lattice sites and the crystal family.
    atol : float, optional
        Absolute tolerance to use for numpy.isclose.  This is used both to
        check for atoms in the equivalent lattice sites and the crystal family.
    check_family : bool, optional
        If True (default), then the crystal family of the unit cell is checked
        to see if the family+setting is one of the 14 standard Bravais space
        lattices.

    Returns
    -------
    bool
        True if the unit cell appears consistent with the indicated Bravais
        lattice setting.  False if atoms are not found at the ideal lattice
        sites, or no Bravais space lattice exists with the family+setting.
    """

    family = ucell.box.identifyfamily(rtol=rtol, atol=atol)

    # Define relative positions based on setting and check crystal family
    if setting == 'p':
        relpos = np.array([[0.0, 0.0, 0.0]])
        families = ['cubic', 'hexagonal', 'tetragonal', 'rhombohedral',
                    'orthorhombic', 'monoclinic', 'triclinic']

    elif setting == 'i':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.5]])
        families = ['orthorhombic', 'tetragonal', 'cubic']

    elif setting == 'f':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.0],
                           [0.5, 0.0, 0.5],
                           [0.0, 0.5, 0.5]])
        families = ['orthorhombic', 'cubic']

    elif setting == 'a':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [0.0, 0.5, 0.5]])
        families = ['monoclinic', 'orthorhombic']

    elif setting == 'b':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.0, 0.5]])
        families = ['monoclinic', 'orthorhombic']

    elif setting == 'c':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [0.5, 0.5, 0.0]])
        families = ['monoclinic', 'orthorhombic']

    elif setting == 't1':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [2.0, 1.0, 1.0],
                           [1.0, 2.0, 2.0]]) / 3.
        families = ['hexagonal']

    elif setting == 't2':
        relpos = np.array([[0.0, 0.0, 0.0],
                           [1.0, 2.0, 1.0],
                           [2.0, 1.0, 2.0]]) / 3.
        families = ['hexagonal']

    else:
        raise ValueError('invalid setting: must be p, i, f, a, b, c, t, t1 or t2')

    # Check that crystal family + setting is a Bravais space lattice
    if check_family and family not in families:
        return False

    pos = ucell.box.position_relative_to_cartesian(relpos)

    atype = None

    for p in pos:
        # Check for atom at pos
        index = index_of_pos(ucell, p, rtol=rtol, atol=atol)
        if np.sum(index) == 0:
            return False
        elif np.sum(index) > 1:
            raise ValueError('Multiple overlapping atoms found')

        # Check that atom's atype is same as other pos
        if atype is None:
            atype = ucell.atoms.atype[index][0]
        elif atype != ucell.atoms.atype[index][0]:
            return False

    return True

def index_of_pos(system, pos, rtol=1e-05, atol=1e-08):
    """
    Searches for atoms in a system to check if they correspond to the given position
    within rounding tolerances and accounting for periodic boundaries.

    Parameters
    ----------
    system : atomman.System
        The system to search.
    pos : array-like object
        The atomic position to search for in the system.
    atol : float, optional
        The absolute tolerance to use with numpy.isclose for identifying matches.
    rtol : float, optional
        The relative tolerance to use with numpy.isclose for identifying matches.

    Returns
    -------
    numpy.ndarray
        Bool values indicating which atoms in system have coordinates at pos.
    """
    # Compute dmag between all atoms in ucell and pos
    dmag = np.array(system.dmag(system.atoms.pos, np.asarray(pos)))

    # Reshape for single atom systems
    if dmag.shape == ():
        dmag = np.array([dmag])

    # Identify which atoms are at pos
    return np.isclose(dmag, 0.0, rtol=rtol, atol=atol)
