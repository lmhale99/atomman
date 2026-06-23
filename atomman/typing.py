from typing import Union

import numpy.typing as npt

from potentials.typing import lammps, unitfloat
from potentials.record.PotentialLAMMPS import PotentialLAMMPS
from potentials.record.PotentialLAMMPSKIM import PotentialLAMMPSKIM

# typing for the LAMMPS potential classes
lammpspotential = Union[PotentialLAMMPS, PotentialLAMMPSKIM]

# typing for lammps executable or lammps library interface
lammpsapp = Union[str, lammps.lammps]

millerindices = Union[str, npt.ArrayLike]

__all__ = ['lammps', 'unitfloat', 'lammpspotential', 'lammpsapp', 'millerindices']