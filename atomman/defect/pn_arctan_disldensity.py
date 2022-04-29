# coding: utf-8
# Standard Python libraries
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# Local imports
from .pn_arctan_disregistry import pn_arctan_disregistry

def pn_arctan_disldensity(x: Optional[npt.ArrayLike] = None,
                          xmax: Optional[float] = None,
                          xstep: Optional[float] = None,
                          xnum: Optional[int] = None,
                          burgers: Union[float, npt.ArrayLike, None] = None,
                          center: float = 0.0,
                          halfwidth: float = 1,
                          normalize: bool = True
                          ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the classic Peierls-Nabarro dislocation density based on an arctan
    disregistry for an array of points x.
        
        ρ(x) = b / π * (ξ / (x^2 + ξ^2)

    Parameters
    ----------
    x : array-like or None, optional
        The exact values of x to use.  Can be specified instead of the xmax,
        xstep, and xnum parameters.
    xmax : float or None, optional
        Maximum value of x to use.  Minimum value is taken as -xmax.  At least
        2 of xmax, xstep, and xnum must be not None.  Default value is None.
    xstep : float or None, optional
        Step size to use between each x value.  At least 2 of xmax, xstep, and
        xnum must be not None.  Default value is None.
    xnum : int or None, optional
        Number of x values to use.  At least 2 of xmax, xstep, and xnum must
        be not None.  Default value is None.
    burgers : float, numpy.ndarray, optional
        The Burgers vector or Burgers vector magnitude for the dislocation. Default
        value is [1, 0, 0].
    center : float
        The x coordinate to center the dislocation at. Default value is 0.0.
    halfwidth : float, optional
        The dislocation halfwidth to use. Default value is 1.
    normalize : bool, optional
        If True (default), the disldensity values will be scaled such that the
        integral over the x range is equal to b.
        
    Returns
    -------
    x : numpy.ndarray
        The x-coordinates for the disregistry values.
    disldensity : numpy.ndarray
        The disldensity vector at each x-coordinate.
    """
    if x is not None:
        if xmax is not None or xstep is not None or xnum is not None:
            raise ValueError('Invalid parameters: x cannot be given with xmax, xstep or xnum.')
    
    else:
        
        # Generate missing x parameters
        if xmax is None:
            if xstep is None or xnum is None:
                raise ValueError('At least two parameters must be given')
            xmax = xstep * (xnum - 1) / 2
        elif xstep is None:
            if xnum is None:
                raise ValueError('At least two parameters must be given')
            xstep = (2 * xmax) / (xnum - 1)
        elif xnum is None:
            xnum = ((2 * xmax) / xstep) + 1

        # Round xnum to int if needed
        if isinstance(xnum, float):
            if np.isclose(xnum, int(xnum + 0.5)):
                xnum = int(xnum + 0.5)
            else:
                raise ValueError('Invalid parameters: xnum or ((2 * xmax) / xstep) not an integer.')

        # Generate x and validate
        x, dx = np.linspace(-1*xmax, xmax, xnum, retstep=True)
        if not np.isclose(dx, xstep):
            raise ValueError('Incompatible parameters: xmax = xstep * (xnum - 1) / 2')
    
    if burgers is None:
        burgers = np.array([1.0, 0.0, 0.0])
    burgers = np.asarray(burgers)
    
    # ρ(x) = b * ξ / (π (x^2 + ξ^2)) 
    disldensity = np.outer(halfwidth / ((x - center)**2 + halfwidth**2), burgers / np.pi)
    
    if normalize is True:
        # Compute the un-normalized disregistry
        disregistry = pn_arctan_disregistry(x=x, burgers=burgers, center=center, halfwidth=halfwidth,
                                            normalize=False)[1]
        # The change in disregistry between xmin and xmax equals the integral of disldensity
        intdisldensity = disregistry[-1] - disregistry[0]
        
        # Rescale disldensity to match the burgers vector
        disldensity = disldensity * np.linalg.norm(burgers) / np.linalg.norm(intdisldensity)
    
    return x, disldensity