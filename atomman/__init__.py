# coding: utf-8

# Standard Python imports
from importlib import resources

# potentials imports
from potentials import build_lammps_potential, settings

# atomman imports
from . import unitconvert
from . import tools
from . import thermo
from . import mep
from . import region
from . import lammps
from .dump import dump, set_dump_styles
from .dump import __all__ as dump_all
from .core import *
from .core import __all__ as core_all
from . import cluster
from . import library
from .library import load_lammps_potential
from .load import load, FileFormatError
from .load import __all__ as load_all
from . import plot
from . import defect

# Set dump styles
set_dump_styles()

# Read version from VERSION file
if hasattr(resources, 'files'):
    __version__ = resources.files('atomman').joinpath('VERSION').read_text(encoding='UTF-8')
else:
    __version__ = resources.read_text('atomman', 'VERSION', encoding='UTF-8').strip()

# Build all list
__all__ = ['__version__', 'load_lammps_potential', 'build_lammps_potential', 'settings',
           'tools', 'thermo', 'mep', 'region', 'lammps', 'dump', 'cluster', 'library',
           'load', 'FileFormatError', 'plot', 'defect']
__all__ += dump_all + core_all + load_all
__all__.sort()

# Define default working units
unitconvert.reset_units(length = 'angstrom', mass = 'amu', energy='eV', charge='e')
