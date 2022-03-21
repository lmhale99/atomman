# coding: utf-8

# Standard Python imports
from typing import Optional
import warnings

# http://www.numpy.org/
import numpy as np

__all__ = ['iscubic', 'ishexagonal', 'istetragonal', 'isrhombohedral',
           'isorthorhombic', 'ismonoclinic', 'istriclinic', 'identifyfamily']

warnmsg = "This is now a method of Box.  The stand-alone function will be depreciated in the next major atomman release."

def iscubic(box, 
            rtol: float = 1e-05,
            atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard cubic cell:
    a = b = c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard cubic cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def ishexagonal(box,
                rtol: float = 1e-05,
                atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard hexagonal cell:
    a = b != c
    alpha = beta = 90
    gamma = 120
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard hexagonal cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 120.0, atol=atol, rtol=rtol))

def istetragonal(box,
                 rtol: float = 1e-05,
                 atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard tetragonal cell:
    a = b != c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard tetragonal cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def isrhombohedral(box,
                   rtol: float = 1e-05,
                   atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard rhombohedral cell:
    a = b = c
    alpha = beta = gamma != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard rhombohedral cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, box.beta, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, box.gamma, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol))

def isorthorhombic(box,
                   rtol: float = 1e-05,
                   atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard orthorhombic cell:
    a != b != c
    alpha = beta = gamma = 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard orthorhombic cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def ismonoclinic(box,
                 rtol: float = 1e-05,
                 atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard monoclinic cell:
    a != b != c
    alpha = gamma = 90
    beta != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard monoclinic cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and np.isclose(box.alpha, 90.0, atol=atol, rtol=rtol)
            and not np.isclose(box.beta, 90.0, atol=atol, rtol=rtol)
            and np.isclose(box.gamma, 90.0, atol=atol, rtol=rtol))

def istriclinic(box,
                rtol: float = 1e-05,
                atol: float = 1e-08) -> bool:
    """
    Tests if a box is consistent with a standard triclinic cell:
    a != b != c
    alpha != 90
    beta != 90
    gamma != 90
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    bool
        True if box is a standard triclinic cell, False otherwise.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    return (not np.isclose(box.a, box.b, atol=atol, rtol=rtol)
            and not np.isclose(box.a, box.c, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, box.beta, atol=atol, rtol=rtol)
            and not np.isclose(box.alpha, box.gamma, atol=atol, rtol=rtol))

def identifyfamily(box, 
                   rtol: float = 1e-05,
                   atol: float = 1e-08) -> Optional[str]:
    """
    Tests if a box is consistent with a standard representation
    of a crystal system cell.
    
    Parameters
    ----------
    box : atomman.Box
        The box object to test.
    rtol : float, optional
        Relative tolerance for testing box parameters. Default value is 1e-5.
    atol : float, optional
        Absolute tolerance for testing box parameters. Default value is 1e-8.
        
    Returns
    -------
    str or None
        'cubic', 'hexagonal', 'tetragonal', 'rhombohedral', 'orthorhombic',
        'monoclinic' or 'triclinic' if it matches any. None if no matches.
        
    Raises
    ------
    ValueError
        If box is not consistent with a standard cell.
    """
    warnings.warn(warnmsg, PendingDeprecationWarning)
    if iscubic(box, rtol=rtol, atol=atol):
        return 'cubic'
    elif ishexagonal(box, rtol=rtol, atol=atol):
        return 'hexagonal'
    elif istetragonal(box, rtol=rtol, atol=atol):
        return 'tetragonal'
    elif isrhombohedral(box, rtol=rtol, atol=atol):
        return 'rhombohedral'
    elif isorthorhombic(box, rtol=rtol, atol=atol):
        return 'orthorhombic'
    elif ismonoclinic(box, rtol=rtol, atol=atol):
        return 'monoclinic'
    elif istriclinic(box, rtol=rtol, atol=atol):
        return 'triclinic'
    else:
        None