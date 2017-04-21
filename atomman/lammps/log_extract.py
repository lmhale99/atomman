import pandas as pd
import warnings
import numpy as np

#atomman imports
from atomman.tools import uber_open_rmode
from DataModelDict import DataModelDict as DM

from .Log import Log


def log_extract(log_info):
    """Parses a LAMMPS screen output/log file and returns a data model containing the information."""
    with warnings.catch_warnings():
        warnings.simplefilter('always')
        warnings.warn('log_extract function is replaced with the Log class', DeprecationWarning)
    
    #Generate a Log object
    log = Log(log_info)
    return log.model()