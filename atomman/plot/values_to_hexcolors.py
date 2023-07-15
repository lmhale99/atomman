# coding: utf-8
from typing import Optional

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

import matplotlib

def values_to_hexcolors(values: npt.ArrayLike,
                        cmap: str,
                        minvalue: Optional[float] = None,
                        maxvalue: Optional[float] = None,
                        reverse: bool = False):
    """
    Takes an array of values and converts them to hex value numbers based on a
    colormap.
    
    Parameters
    ----------
    values : Array-like object
        The values to use as the basis for coloring.
    cmap : str or matplotlib colormap
        The matplotlib colormap to use for identifying the colors.
    minvalue : float or None, optional
        The minimum limit value to use for the color scale.  If None (default)
        then the minimum value in values will be used.
    maxvalue : float or None, optional
        The maximum limit value to use for the color scale.  If None (default)
        then the maximum value in values will be used.
    reverse : bool, optional
        If set to True the colormap will be reversed.
    """
    # Get colormap
    if isinstance(cmap, str):
        cmap = matplotlib.colormaps[cmap]

    # Set min and max values
    if minvalue is None:
        minvalue = values.min()
    if maxvalue is None:
        maxvalue = values.max()

    # Scale values to range from 0 to 1
    if minvalue == maxvalue:
        scalevalues = np.full_like(values, 0.5)
    else:
        scalevalues = (values - minvalue) / (maxvalue - minvalue)

    # Reverse values
    if reverse:
        scalevalues = 1 - scalevalues

    # Define converter from value to hexcolor and vectorize
    def hexcolors(value):
        return matplotlib.colors.rgb2hex(cmap(value))
    vhexcolors = np.vectorize(hexcolors)

    # Convert scaled values to hexcolors
    colors = vhexcolors(scalevalues)

    return colors
