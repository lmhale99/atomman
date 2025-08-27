# coding: utf-8

class LammpsError(Exception):
    """
    Class for passing LAMMPS errors to Python
    """
    def __init__(self, *args, **kwargs):
        
        message = args[0]
        try:
            # Crop message to only what appears after ERROR(something):
            error_index = message.rindex('ERROR')
            crop_index = message[error_index:].index(':') + error_index + 1
            message = message[crop_index:].strip()
        except:
            pass
        args = (message, ) + args[1:]

        super().__init__(*args, **kwargs)

from . import style
from .Log import Log
from .NEBLog import NEBLog
from .Potential import Potential, PotentialLAMMPS, PotentialLAMMPSKIM
from .run import run
from .normalize import normalize
from .checkversion import checkversion
from .seed import seed, seedmax, newseed
__all__ = ['LammpsError', 'style', 'run', 'normalize', 'Potential', 'Log',
           'PotentialLAMMPS', 'PotentialLAMMPSKIM', 'NEBLog', 'checkversion',
           'seed', 'seedmax', 'newseed']