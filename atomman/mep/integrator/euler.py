# coding: utf-8

def euler(ratefxn, coord, timestep, **kwargs):
    """
    Performs Euler ODE integration for a timestep.
    
    Parameters
    ----------
    ratefxn : function
        The rate function to use.
    coord : array-like object
        The coordinate(s) of the last timestep.
    timestep : float
        The timestep value to use.
    
    Returns
    -------
    array-like object
        The coordinate(s) moved forward by timestep.
    """
    
    return coord + timestep * ratefxn(coord, **kwargs)