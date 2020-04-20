# coding: utf-8
# Standard Python libraries
from importlib import resources

# Read version from VERSION file
__version__ = resources.read_text('atomman', 'VERSION').strip()

# atomman imports
from . import unitconvert
from . import tools
from . import region
from . import lammps
from .dump import *
from .dump import __all__ as dump_all
from .core import *
from .core import __all__ as core_all
from .load import *
from .load import __all__ as load_all
from . import plot
from . import defect


__all__ = ['__version__'] + dump_all + core_all + load_all
__all__.sort()

# Define default working units
unitconvert.reset_units(length = 'angstrom', mass = 'amu', energy='eV', charge='e')