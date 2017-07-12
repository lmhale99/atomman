# Define rootdir
import os
rootdir = os.path.dirname(os.path.abspath(__file__))

# Read version from VERSION file
with open(os.path.join(rootdir, 'VERSION')) as version_file:
    __version__ = version_file.read().strip()

from .core import *
from . import unitconvert
from . import defect
from . import lammps
from . import tools
from . import convert
from . import plot

unitconvert.reset_units(length = 'angstrom', mass = 'amu', energy='eV', charge='e')
