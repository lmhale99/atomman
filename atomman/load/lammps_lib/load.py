# Standard Python libraries
from typing import Union

import numpy as np

# atomman imports
import atomman.unitconvert as uc
from ... import Atoms, Box, System
from ...lammps import style
from ...typing import lammps

def load(lmp: lammps.lammps,
         symbols: Union[str, list, None] = None,
         lammps_units: str = 'metal',) -> System:
    """
    Extract system information from a LAMMPS library object and build
    an atomman.System

    NOTE: This currently only reads in box information, atom types and
    atom positions!  Support for additional per-atom properties to come...

    Parameters
    ----------
    lmp : lammps.lammps
        The LAMMPS object from which the system information will be retrieved.
    symbols : str or list, optional
        Allows the list of element symbols to be assigned during loading.
    lammps_units : str
        The LAMMPS units option associated with the parameters.  Default value
        is 'metal'.
    """
    lammps_unit = style.unit(lammps_units)

    # Get box parameters and build Box
    lo, hi, xy, yz, xz, periodicity, box_change = lmp.extract_box()
    xlo = uc.set_in_units(lo[0], lammps_unit['length'])
    xhi = uc.set_in_units(hi[0], lammps_unit['length'])
    ylo = uc.set_in_units(lo[1], lammps_unit['length'])
    yhi = uc.set_in_units(hi[1], lammps_unit['length'])
    zlo = uc.set_in_units(lo[2], lammps_unit['length'])
    zhi = uc.set_in_units(hi[2], lammps_unit['length'])
    xy = uc.set_in_units(xy, lammps_unit['length'])
    xz = uc.set_in_units(xz, lammps_unit['length'])
    yz = uc.set_in_units(yz, lammps_unit['length'])
    box = Box(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi,
              xy=xy, yz=yz, xz=xz)
    pbc = np.array(periodicity, dtype=bool)

    # Get atom info and build Atoms
    natoms = lmp.get_natoms()
    atype = lmp.numpy.extract_atom('type', nelem=natoms, dim=1)
    pos = lmp.numpy.extract_atom('x', nelem=natoms, dim=3)
    pos = uc.set_in_units(pos, lammps_unit['length'])
    atoms = Atoms(natoms=natoms, atype=atype, pos=pos, safecopy=True)

    # Build system
    system = System(atoms=atoms, box=box, pbc=pbc, symbols=symbols)

    return system