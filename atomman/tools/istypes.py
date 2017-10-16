# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

# atomman imports
from ..compatibility import long

__all__ = ['is_int', 'is_dtype_int', 'is_bool', 'is_dtype_bool']

def is_int(value):
    """
    Determines if a number is of any integer type (standard or numpy based).
    
    Parameters
    ----------
    value : any
        The value to test
        
    Returns
    -------
    bool
        True if value is int, long, or 0D numpy array with an int dtype.  
        False otherwise (including booleans).
    """
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, long)): 
        return True
    elif isinstance(value, np.ndarray) and value.shape == ():
        return is_dtype_int(value.dtype)
    else:
        return False
        
def is_dtype_int(dtype):
    """
    Determines if a numpy array's dtype is an integer type (excluding
    booleans).
    
    Parameters
    ----------
    dtype: numpy.dtype
        The datatype to test.
    
    Returns
    -------
    bool
        True if dtype is any int type.
        False otherwise (including dtype == bool).
    """
    if dtype == bool:
        return False
    else:
        try:
            if np.issubdtype(dtype.type, int):
                return True
            else:
                return False
        except:
            return False

def is_bool(value):
    """
    Determines if a number is of any boolean type (standard or numpy based).
    
    Parameters
    ----------
    value : any
        The value to test
        
    Returns
    -------
    bool
        True if value is bool, or 0D numpy array with an bool dtype.
        False otherwise (including booleans).
    """
    if isinstance(value, bool): 
        return True
    elif isinstance(value, np.ndarray) and value.shape == ():
        return is_dtype_bool(value.dtype)
    else:
        return False
          
def is_dtype_bool(dtype):
    """
    Determines if the dtype of a numpy array is of any bool type.
    
    Parameters
    ----------
    dtype: numpy.dtype
        The datatype to test.
    
    Returns
    -------
    bool
        True if dtype is bool.
        False otherwise.
    """
    return dtype == bool