# coding: utf-8

# Standard Python libraries
import ast
from typing import Optional

# https://pypi.python.org/pypi/numericalunits
import numericalunits as nu

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# yabadaba imports
from yabadaba import unitconvert as unitconvert_obj

# NOTES: This module now serves as a wrapper around yabadaba.unitconvert!


unit = unitconvert_obj.unit

def build_unit():
    """
    Saves numericalunits attributes to global dictionary unit so the values
    can be retrieved by their string names.
    """
    unitconvert_obj.build_unit()
    #unit.update(unitconvert_obj.unit)

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
    #unit.update(unitconvert_obj.unit)

def set_literal(term: str) -> npt.ArrayLike:
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

def set_in_units(value: npt.ArrayLike,
                 units: str) -> npt.ArrayLike:
    """
    Convert value from specified units to working units.
    
    Parameters
    ----------
    value : array-like object
        A numerical value or list/array of values.
    units : str
        The units that value is in.
        
    Returns
    -------
    float or numpy.ndarray
        The given value converted from the specified units to working units.
    """
    return unitconvert_obj.set_in_units(value=value, units=units)

def get_in_units(value: npt.ArrayLike,
                 units: str) -> npt.ArrayLike:
    """
    Convert value from working units to specified units.
    
    Parameters
    ----------
    value : array-like object
        A numerical value or list/array of values.
    units : str
        The units to convert value to (from working units).
        
    Returns
    -------
    float or numpy.ndarray
        The given value converted to the specified units from working units.
    """
    return unitconvert_obj.get_in_units(value=value, units=units)

def value_unit(term: dict) -> npt.ArrayLike:
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
    
def error_unit(term: dict) -> npt.ArrayLike:
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
    
def model(value: npt.ArrayLike,
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
        values to interpret.
        
    Returns
    -------
    float
        The scaling factor for converting numbers in the given units to
        working units. If units is None or == 'scaled', then this value is
        1.0.
    """
    return unitconvert_obj.parse(units=units)