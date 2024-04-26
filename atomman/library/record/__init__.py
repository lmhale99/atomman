# coding: utf-8

from potentials.record import Record, load_record, recordmanager
__all__ = ['Record', 'load_record', 'recordmanager']
__all__.sort()

# Add the modular Record styles
recordmanager.import_style('crystal_prototype', '.CrystalPrototype', __name__)
recordmanager.import_style('relaxed_crystal', '.RelaxedCrystal', __name__)
recordmanager.import_style('reference_crystal', '.ReferenceCrystal', __name__)
recordmanager.import_style('free_surface', '.FreeSurface', __name__)
recordmanager.import_style('stacking_fault', '.StackingFault', __name__)
recordmanager.import_style('point_defect', '.PointDefect', __name__)
recordmanager.import_style('dislocation', '.Dislocation', __name__)
recordmanager.import_style('grain_boundary', '.GrainBoundary', __name__)

