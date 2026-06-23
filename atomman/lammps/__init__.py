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
from .Potential import Potential
from .run import run, run_libtest, restart_check, read_logs
from .normalize import normalize
from .checkversion import checkversion, versiondate
from .seed import seed, seedmax, newseed
from .LAMMPS import LAMMPS, LAMMPSEXE, LAMMPSLIB, LAMMPSobj

__all__ = ['LammpsError', 'style', 'run', 'run_libtest', 'normalize', 'Potential', 'Log',
           'NEBLog', 'checkversion', 'seed', 'seedmax', 'newseed', 'versiondate',
           'restart_check', 'read_logs', 'LAMMPS', 'LAMMPSEXE', 'LAMMPSLIB',
           'LAMMPSobj', 'anylammps']