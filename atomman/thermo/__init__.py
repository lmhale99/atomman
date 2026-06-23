from .RDF import RDF
from .StructureFactor import StructureFactor
from .PartialStructureFactor import PartialStructureFactor
from .IdealGas import IdealGas
from .EinsteinSolid import EinsteinSolid
from .UhlenbeckFordModel import UhlenbeckFordModel

__all__ = ['RDF', 'StructureFactor', 'PartialStructureFactor', 'IdealGas', 'EinsteinSolid', 'UhlenbeckFordModel']

from potentials.record import recordmanager

# Add the modular Record styles
recordmanager.import_style('rdf', '.RDF', __name__)
recordmanager.import_style('structure_factor', '.StructureFactor', __name__)
recordmanager.import_style('partial_structure_factor', '.PartialStructureFactor', __name__)