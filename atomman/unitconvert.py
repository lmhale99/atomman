# coding: utf-8

# Standard Python libraries
from typing import Optional, Union

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# yabadaba imports
from yabadaba import unitconvert as unitconvert_obj

# NOTES: This module now serves as a wrapper around yabadaba.unitconvert!

# Local imports
from . import typing as amt


unit = unitconvert_obj.unit

def build_unit():
    """
    Saves numericalunits attributes to global dictionary unit so the values
    can be retrieved by their string names.
    """
    unitconvert_obj.build_unit()

def reset_units(seed: Optional[int] = None, **kwargs):
    """
    Extends numericalunits.reset_units() by allowing for working units to be
    defined.  If no working units are specified, then random working units are
    used just like the default numericalunits behavior.  Otherwise, use the
    specified working units and SI.
    
    Parameters
    ----------
    seed : int, optional
        random number seed to use in generating random working units.
        seed='SI' will use SI units.  Cannot be given with the other
        parameters.
    length : str, optional
        Unit of length to use for the working units.
    mass : str, optional
        Unit of mass to use for the working units.
    time : str, optional
        Unit of time to use for the working units.
    energy : str, optional
        Unit of energy to use for the working units.
    charge : str, optional
        Unit of charge to use for the working units.
        
    Raises
    ------
    ValueError
        If seed is given with any other parameters, or if more than four of
        the working unit parameters are given.
    """
    unitconvert_obj.reset_units(seed=seed, **kwargs)

def set_literal(term: str) -> np.ndarray:
    """
    Convert string 'value unit' to numbers in working units.
    
    Parameters
    ----------
    term : str
        String containing value and associated unit. If unit is not given,
        then the value is converted to a float and assumed to be in working
        units.
        
    Returns
    -------
    float or numpy.ndarray
        The numerical value of term in working units.
        
    Raises
    ------
    ValueError
        If no valid float value can be parsed.
    """
    return unitconvert_obj.set_literal(term)

def set_in_units(value: amt.unitfloat,
                 units: Optional[str] = None) -> np.ndarray:
    """
    Parameters
    ----------
    value : str, float or array-like object
        A numerical value or list/array of values.  String values will be
        interpreted using set_literal() and may therefore optionally contain
        units information.  If a string value contains units information, do
        not give it separately to the units parameter!
    units : str or None
        The units that the value is in.  The default value of None will either
        do no conversion or use units defined in a string value.
        
    Returns
    -------
    float or numpy.ndarray
        The given value converted from the specified units to working units.
    """
    return unitconvert_obj.set_in_units(value=value, units=units)

def get_in_units(value: Union[float, npt.ArrayLike],
                 units: Optional[str] = None) -> np.ndarray:
    """
    Convert value from working units to specified units.
    
    Parameters
    ----------
    value : float or array-like object
        A numerical value or list/array of values.
    units : str or None
        The units to convert value to (from working units).  A value
        of None will do no conversion.
        
    Returns
    -------
    float or numpy.ndarray
        The given value converted to the specified units from working units.
    """
    return unitconvert_obj.get_in_units(value=value, units=units)

def value_unit(term: dict) -> np.ndarray:
    """
    Reads numerical value from dictionary containing 'value' and 'unit' keys.
    
    Parameters
    ----------
    term : dict
        Dictionary containing 'value' and 'unit' keys.
        
    Returns
    -------
    float or numpy.ndarray
        The result of calling set_in_units() by passing the dictionary keys 
        'value' and 'unit' as parameters.
    
    """
    return unitconvert_obj.value_unit(term=term)
    
def error_unit(term: dict) -> np.ndarray:
    """
    Reads numerical error from dictionary containing 'error' and 'unit' keys.
    
    Parameters
    ----------
    term : dict
        Dictionary containing 'error' and 'unit' keys.
        
    Returns
    -------
    float or numpy.ndarray
        The result of calling set_in_units() by passing the dictionary keys 
        'error' and 'unit' as parameters.
    
    """
    return unitconvert_obj.error_unit(term=term)
    
def model(value: Union[float, npt.ArrayLike],
          units: Optional[str] = None,
          error: Optional[npt.ArrayLike] = None) -> DM:
    """
    Generates DataModelDict representation of data.
    
    Parameters
    ----------
    value : array-like object
        A numerical value or list/array of values.
    units : str, optional
        The units to convert value to (from working units).
    error : array-like object or None, optional
        A value error to include.  If given, must be the same
        size/shape as value.
    
    Returns
    -------
    DataModelDict
        Model representation of the value(s).
    """
    return unitconvert_obj.model(value=value, units=units, error=error)

def parse(units: Optional[str]) -> float:
    """
    Convert units as strings (or None) into scaling numbers.  This function
    allows for complex unit definitions with operators:
    
    - '()' for defining order of operations
    - '*' for multiplication.
    - '/' for division.
    - '^' for powers.
    
    Parameters
    ----------
    units : str or None
        String consisting of defined unit names, operators, and numerical 
        values to interpret.  If units is None or == 'scaled', then the returned
            value will be 1.0, i.e. no unit conversion.
        
    Returns
    -------
    float
        The scaling factor for converting numbers in the given units to
        working units.
    """
    return unitconvert_obj.parse(units=units)