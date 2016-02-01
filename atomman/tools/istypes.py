def isint(value):
    """Determines if a number is of any integer type (standard or numpy based)."""
    if isinstance(value, (int, long)): 
        return True
    elif isinstance(value, np.ndarray) and len(value) == 0:
        return is_dtype_int(value.dtype)
    else:
        return False
        
def is_dtype_int(dtype):
    """Determines if the dtype of a numpy array is of any integer type (standard or numpy based)."""
    if dtype == int or dtype == long: 
        return True
    else:
        try:
            if issubclass(dtype.type, np.int):
                return True
            else:
                return False
        except:
            return False

def isbool(value):
    """Determines if a number is of any boolean type (standard or numpy based)."""
    if isinstance(value, bool): 
        return True
    elif isinstance(value, np.ndarray) and len(value) == 0:
        return is_dtype_bool(value.dtype)
    else:
        return False
          
def is_dtype_bool(dtype):
    """Determines if the dtype of a numpy array is of any bool type (standard or numpy based)."""
    if dtype == bool: 
        return True
    else:
        try:
            if issubclass(dtype.type, np.bool):
                return True
            else:
                return False
        except:
            return False