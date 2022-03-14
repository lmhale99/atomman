# coding: utf-8

# Standard Python imports
from typing import Callable, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

def euler(ratefxn: Callable,
          coord: npt.ArrayLike,
          timestep: float,
          **kwargs) -> Union[float, np.ndarray]:
    """
    Performs Euler ODE integration for a timestep.
    
    Parameters
    ----------
    ratefxn : function
        The rate function to use.  Should be a function of coord.
    coord : array-like object
        The coordinate(s) of the last timestep.
    timestep : float
        The timestep value to use.
    **kwargs : any
        Any extra keyword parameters to pass on to ratefxn.
    
    Returns
    -------
    float or numpy.ndarray
        The coordinate(s) moved forward by timestep.
    """
    
    return coord + timestep * ratefxn(coord, **kwargs)