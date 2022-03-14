# coding: utf-8

# Standard Python imports
from typing import Callable, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

def rungekutta(ratefxn: Callable,
               coord: npt.ArrayLike,
               timestep: float,
               **kwargs) -> Union[float, np.ndarray]:
    """
    Performs Runge-Kutta ODE integration for a timestep.
    
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
    
    k1 = timestep * ratefxn(coord, **kwargs)
    k2 = timestep * ratefxn(coord - 0.5 * k1, **kwargs)
    k3 = timestep * ratefxn(coord - 0.5 * k2, **kwargs)
    k4 = timestep * ratefxn(coord - k3, **kwargs)
    
    return coord + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6