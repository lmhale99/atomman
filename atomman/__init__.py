# coding: utf-8

# Read version from VERSION file
try:
    from importlib import resources
except:
    from pathlib import Path
    with open(Path(Path(__file__).resolve().parent, 'VERSION')) as version_file:
        __version__ = version_file.read().strip()
else:
    __version__ = resources.read_text('atomman', 'VERSION').strip()

# potentials imports
from potentials import build_lammps_potential, Settings

# atomman imports
from . import unitconvert
from . import tools
from . import mep
from . import region
from . import lammps
from .dump import *
from .dump import __all__ as dump_all
from .core import *
from .core import __all__ as core_all
from . import cluster
from . import library
from .library import load_lammps_potential
from .load import *
from .load import __all__ as load_all
from . import plot
from . import defect


__all__ = ['__version__'] + dump_all + core_all + load_all
__all__ += ['load_lammps_potential', 'build_lammps_potential']
__all__.sort()

# Define default working units
unitconvert.reset_units(length = 'angstrom', mass = 'amu', energy='eV', charge='e')