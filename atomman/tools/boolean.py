# coding: utf-8

# Standard Python libraries
from typing import Union

def boolean(value: Union[str, bool]) -> bool:
    """
    Allows conversion of strings to Booleans.
    
    Parameters
    ----------
    value : str or bool
        If str, then 'true' and 't' become True and 'false' and 'f' become
        false. If bool, simply return the value.
        
    Returns
    -------
    bool
        Equivalent bool of value.
        
    Raises
    ------
    ValueError
        If value is unrecognized.
    """
    
    # Pass Boolean values through without changing
    if value is True:
        return True
    elif value is False:
        return False
    
    # Convert strings
    elif value.lower() in ['true', 't']:
        return True
    elif value.lower() in ['false', 'f']:
        return False
    
    # Issue error for invalid string
    else:
        raise ValueError('Invalid Boolean string')