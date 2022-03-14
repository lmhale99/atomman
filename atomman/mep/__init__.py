# coding: utf-8

# Standard Python imports
from typing import Callable, Optional

# http://www.numpy.org/
import numpy.typing as npt

# atomman imports
from . import gradient, integrator
from .BasePath import BasePath
from .ISMPath import ISMPath

__all__ = ['gradient', 'integrator', 'BasePath', 'ISMPath', 'create_path']

def create_path(coord: npt.ArrayLike,
                energyfxn: Callable,
                style: str = 'ISM',
                gradientfxn: str = 'cdiff',
                gradientkwargs: Optional[dict] = None,
                integratorfxn: str = 'rk') -> BasePath:
    """
    Generic function call for initializing a Path object.
    
    Parameters
    ----------
    coord : array-like object
        The list of coordinates associated with the points along the path.
    energyfxn : function
        The function that evaluates the energy associated with the different
        point coordinates.
    style : str, optional
        The relaxation style to associate with the path.  Determines the
        subclass to build.  Default value is "ISM" for the improved
        string method.
    gradientfxn : function, optional
        The function to use to estimate the gradient of the energy.  Default
        value of 'cdiff' will use atomman.mep.gradient.central_difference
    gradientkwargs : dict, optional
        The keyword arguments (i.e. settings) to use with the gradientfxn.
        Default is an empty dictionary, i.e. default settings of gradientfxn.
    integratorfxn : str or function, optional
        The function to use to integrate relaxation steps.  Default value of
        'rk' will use atomman.mep.integrator.rungekutta.

    Returns
    -------
    BasePath
        A Path object.

    Raises
    ------
    ValueError
        If the given style is invalid.
    """
    if style == 'ISM' or style == 'improved_string_method':
        return ISMPath(coord, energyfxn, gradientfxn=gradientfxn,
                       gradientkwargs=gradientkwargs,
                       integratorfxn=integratorfxn)
    else:
        raise ValueError('Unknown style')

