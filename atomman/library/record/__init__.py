import sys

from potentials.record import Record, load_record, recordmanager
__all__ = ['Record', 'load_record', 'recordmanager']

#### Full record styles - include in recordmanager ####

# Import CrystalPrototype
try:
    from .CrystalPrototype import CrystalPrototype
except Exception as e:
    recordmanager.failed_styles['crystal_prototype'] = '%s: %s' % sys.exc_info()[:2]
else:
    recordmanager.loaded_styles['crystal_prototype'] = CrystalPrototype
    __all__.append('CrystalPrototype')

# Import RelaxedCrystal
try:
    from .RelaxedCrystal import RelaxedCrystal
except Exception as e:
    recordmanager.failed_styles['relaxed_crystal'] = '%s: %s' % sys.exc_info()[:2]
else:
    recordmanager.loaded_styles['relaxed_crystal'] = RelaxedCrystal
    __all__.append('RelaxedCrystal')

__all__.sort()