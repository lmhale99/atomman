# coding: utf-8

# Standard Python libraries
from warnings import warn
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://www.scipy.org/
from scipy.interpolate import griddata

# https://matplotlib.org/
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm

# atomman imports
import atomman.unitconvert as uc
from . import get_prop_values
from .. import Box, System
from ..tools import axes_check

__all__ = ['interpolate_contour']

def interpolate_contour(system: System,
                        prop_name: str,
                        prop_index: Union[int, tuple, None] = None,
                        prop_magnitude: bool = False,
                        prop: Optional[npt.ArrayLike] = None,  
                        plotxaxis: str = 'x',
                        plotyaxis: str = 'y',
                        xlim: Optional[tuple] = None,
                        ylim: Optional[tuple] = None,
                        zlim: Optional[tuple] = None,
                        vlim: Optional[tuple] = None,
                        xbins: int = 200,
                        ybins: int = 200,
                        dots: bool = True,
                        colorbar: bool = True,
                        vzero: bool = True,
                        title: Union[bool, str] = True,
                        length_unit: str = 'angstrom',
                        prop_unit: Optional[str] = None,
                        cmap: str = 'jet',
                        fill_value: float = np.nan,
                        figsize: Union[int, float, tuple, None] = None,
                        matplotlib_axes: Optional[plt.axes] = None,
                       ) -> Tuple[float, float, Optional[plt.figure]]:
    """
    Creates a contour plot of a system's per-atom properties by interpolating
    properties between atoms.

    Parameters
    ----------
    system : atomman.System
        The system with the per-atom property that you want to plot.
    propname : str
        The name of the per-atom property that you want to plot.
    prop : array-like object, optional
        Values for the per-atom property to plot.  If not given, values will
        be taken as the "name" property of system.
    propindex : int or tuple, optional
        Specifies which component of a multidimensional property to plot.  Not
        needed if the property is scalar.
    magnitude : bool, optional
        If True, plots the per-atom magnitude of a vector property.  Cannot be
        combined with index.  Default value is False.
    plotxaxis : str or array-like object, optional
        Indicates the Cartesian direction associated with the system's atomic
        coordinates to align with the plotting x-axis.  Values are either 3D
        unit vectors, or strings 'x', 'y', or 'z' for the Cartesian axes
        directions.  plotxaxis and plotyaxis must be orthogonal.  Default value
        is 'x' = [1, 0, 0].
    plotyaxis : str or array-like object, optional
        Indicates the Cartesian direction associated with the system's atomic
        coordinates to align with the plotting y-axis.  Values are either 3D
        unit vectors, or strings 'x', 'y', or 'z' for the Cartesian axes
        directions.  plotxaxis and plotyaxis must be orthogonal.  Default value
        is 'y' = [0, 1, 0].
    xlim : tuple, optional
        The minimum and maximum coordinates along the plotting x-axis to
        include in the fit.  Values are taken in the specified length_unit.
        If not given, then the limits are set based on min and max atomic
        coordinates along the plotting axis.
    ylim : tuple, optional
        The minimum and maximum coordinates along the plotting y-axis to
        include in the fit.  Values are taken in the specified length_unit.
        If not given, then the limits are set based on min and max atomic
        coordinates along the plotting axis.
    zlim : tuple, optional
        The minimum and maximum coordinates normal to the plotting axes
        (i.e. plotxaxis X plotyaxis) to include in the fit.  Values are taken
        in the specified length_unit.  If not given, then the limits are set
        based on min and max atomic coordinates along the axis.
    vlim : tuple, optional
        Allows for the color range to be specified.  If not give, the range
        will be auto-selected based on the vzero parameter. 
    xbins : int, optional
        Specifies the number of interpolation bins to use along the plotting 
        x-axis.  Default value is 200.
    ybins : int, optional
        Specifies the number of interpolation bins to use along the plotting
        y-axis.  Default value is 200.
    dots : bool, optional
        If True, then the positions of the atoms are shown as circles.
        Default value is True.
    vzero : bool, optional
        Specifies how to auto-select the vlim values if vlim is not directly
        given. vzero=True (default) will center the range around 0, i.e.
        vlim[0] = -vlim[1].  vzero=False will select vlim[0] and vlim[1] 
        independently based on the the property values.
    colorbar: bool, optional
        If True (default) a colorbar will be added to the plot.
    title: bool or str, optional
        If True (default), a plot title will be added matching the property
        name plus index/magnitude info.  If False, no title will be used.
        If a string is given, then that will be used instead of the property
        name.
    length_unit : str, optional
        The units of length to use for the plotting x- and y-axes.  Default
        value is 'angstrom'.
    prop_unit : str or None, optional
        The units to use for the property value being plotted.  Default value
        is None, in which no unit conversion is applied.
    cmap : str, optional
        The name of the matplotlib colormap to use.  Default value is 'jet'.
    fill_value: float, optional
        Value used to fill in for grid points failed to interpolate in the fit.
        If not given, then the default is np.nan, which may cause an error in
        plotting for too narrow xlim and ylim settings.
    figsize : int, float or tuple, optional
        Specifies the size of the figure to create in inches.  Note the plot itself
        is always made based on an equal aspect ratio of the coordinates.  If a
        single value is given, then the size of the larger dimension will be set to
        this and the smaller dimension will be scaled accordingly.  The width or
        height can be directly set to a specific value by giving (w, None) or
        (None, w), respectively.  If you give values for both width and height
        the width value affects the size of the plot and height affects the size
        of the colorbar if one is included.
    matplotlib_axes : matplotlib.Axes.axes, optional
        An existing matplotlib axes object. If given, the differential displacement
        plot will be added to the specified axes of an existing figure.  This
        allows for subplots to be constructed.  Note that figsize will be ignored
        as the figure would have to be created beforehand and no automatic
        optimum scaling of the figure's dimensions will occur.
    
    Returns
    -------
    intsum : float
        The area integrated sum of the property over the plotted region.
    avsum : float
        The average property value taken across all plotting bins.
    figure : matplotlib.pyplot.figure
        The generated figure object is returned if matplotlib_axes is not given.
    """

    # Interpret plot axis values
    plotxaxis = __plotaxisoptions(plotxaxis)
    plotyaxis = __plotaxisoptions(plotyaxis)

    # Build transformation matrix, T, from plot axes.
    T = axes_check([plotxaxis, plotyaxis, np.cross(plotxaxis, plotyaxis)])

    # Extract positions and transform using T
    pos = uc.get_in_units(np.inner(system.atoms.pos, T), length_unit)

    # Set plot limits
    if xlim is None:
        xlim = (pos[:, 0].min(), pos[:, 0].max())
    if ylim is None:
        ylim = (pos[:, 1].min(), pos[:, 1].max())
    if zlim is None:
        zlim = (pos[:, 2].min(), pos[:, 2].max())

    # Define box for identifying only points inside
    plotbox = Box(xlo = xlim[0] - 5, xhi = xlim[1] + 5,
                  ylo = ylim[0] - 5, yhi = ylim[1] + 5,
                  zlo = zlim[0], zhi = zlim[1])
    in_bounds = plotbox.inside(pos)

    # Manage property values
    prop_name, prop = get_prop_values(system, prop_name, prop_index,
                                      prop_magnitude, prop, prop_unit)

    # Extract plotting coordinates and values
    x = pos[in_bounds, 0]
    y = pos[in_bounds, 1]
    c = prop[in_bounds]

    # Generate plotting grid
    grid = __grid_interpolate(x, y, c, xbins, ybins, xlim, ylim, fill_value)

    # Compute intsum and avsum values
    intsum = np.sum(grid)
    avsum = intsum / (xbins-1) / (ybins-1)
    intsum = intsum * (xlim[1] - xlim[0]) / (xbins-1) * (ylim[1] - ylim[0]) / (ybins-1)
    returning = [intsum, avsum]

    lx = xlim[1] - xlim[0]
    ly = ylim[1] - ylim[0]

    # Create figure if needed
    if matplotlib_axes is None:
        if figsize is None:
            figsize = [8, None]

        if isinstance(figsize, (int, float)):

            if lx > ly:
                width = figsize
                height = width * ly / lx
            else:
                height = figsize
                width = height * lx / ly
            figsize = (width, height)

        else:
            figsize = list(figsize)
            assert len(figsize) == 2, 'figsize must be a single value or list/tuple of two values'
            if figsize[0] is None:
                assert figsize[1] is not None, 'both figsize terms cannot be None'
                figsize[0] = figsize[1] * lx / ly
            elif figsize[1] is None:
                figsize[1] = figsize[0] * ly / lx

        fig = plt.figure(figsize=figsize)
        matplotlib_axes = fig.add_subplot(111)
        returning.append(fig)
    else:
        fig = matplotlib_axes.get_figure()

    # Get cmap from str name
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)

    if vlim is None:
        vlim = __vlim(grid, vzero)

    # Generate plot
    im = matplotlib_axes.imshow(grid, extent = [xlim[0], xlim[1], ylim[1], ylim[0]], 
                                cmap = cmap, norm = plt.Normalize(vmin=vlim[0], vmax=vlim[1]))

    # Add pretty colorbar
    if colorbar:
        ticks = np.linspace(vlim[0], vlim[1], 11, endpoint=True)
        cbar = fig.colorbar(im, ax=matplotlib_axes, fraction=0.05, pad=0.05, ticks=ticks)
        cbar.ax.tick_params(labelsize=15)

    # Make axes pretty
    matplotlib_axes.axis([xlim[0], xlim[-1], ylim[0], ylim[-1]])
    #matplotlib_axes.tick_params(labelsize=15)

    # Add title
    if isinstance(title, str):
        prop_name = title
        title = True
    if title is True:
        matplotlib_axes.set_title(prop_name, size=30)

    # Add dots
    if dots:
        __adddots(matplotlib_axes, x, y, xlim, ylim)

    return returning

def __plotaxisoptions(plotaxis: Union[str, npt.ArrayLike]) -> np.ndarray:
    """Converts str plotaxis options to arrays as needed"""

    # Give numeric values for str plot axis terms
    if plotaxis == 'x':
        plotaxis = [1.0, 0.0, 0.0]
    elif plotaxis == 'y':
        plotaxis = [0.0, 1.0, 0.0]
    elif plotaxis == 'z':
        plotaxis = [0.0, 0.0, 1.0]

    # Convert to numpy array
    return np.asarray(plotaxis, dtype=float)

def __grid_interpolate(x: npt.ArrayLike,
                       y: npt.ArrayLike,
                       v: npt.ArrayLike,
                       xbins: int = 50,
                       ybins: int = 50,
                       xlim: Optional[tuple] = None,
                       ylim: Optional[tuple] = None,
                       fill_value: float = np.nan) -> np.ndarray:
    """
    Generates 2D grid of property values by interpolating between measured
    values.

    Parameters
    ----------
    x : array-like object
        The list of x-coordinate values.
    y : array-like object
        The list of y-coordinate values.
    v : array-like object
        The list of property values associated with each (x, y) set.
    xbins : int, optional
        The number of bins to use in interpolating between the x coordinates.
        Default value is 50.
    ybins : int, optional
        The number of bins to use in interpolating between the y coordinates.
        Default value is 50.
    xlim : tuple, optional
        The minimum and maximum coordinates along the plotting x-axis to
        include in the fit.
    ylim : tuple, optional
        The minimum and maximum coordinates along the plotting y-axis to
        include in the fit.
    fill_value: float, optional
        Value used to fill in for grid points failed to interpolate in the fit.
        If not given, then the default is np.nan.

    Returns
    -------
    grid : numpy.ndarray
        The interpolated values in a grid map array.
    """
    # Generate 1D interpolation points
    xi = np.linspace(xlim[0], xlim[1], num=xbins)
    yi = np.linspace(ylim[0], ylim[1], num=ybins)

    # Generate 2D grid points
    x0, y0 = np.meshgrid(xi, yi)

    # Interpolate values to grid points
    grid = griddata((x, y), v, (x0, y0), fill_value=fill_value)
    if np.any(np.isnan(grid)):
        warn("Given xlim and ylim are too broad to interpolate. Consider shrinking xlim and ylim or set an appropriate value in fill_value.")

    return grid

def __vlim(grid: np.ndarray,
           vzero: bool = True,
           scale: float = 1) -> Tuple[float, float]:
    """
    Identifies limits to use for the property values
    """
    # Set vmin = -vmax if vzero is True
    if vzero is True:
        vmax = abs(grid).max()
        if vmax != 0.0:
            rounder = np.floor(np.log10(vmax))
            vmax = np.around(2 * vmax / 10.**rounder) * 10.**rounder / 2.
        else:
            vmax = 1e-15
        vmin = -vmax

    # Set vmin, vmax if vzero is False
    else:
        vmax = grid.max()
        vmin = grid.min()

        if abs(grid).max() != 0.0:
            rounder = np.floor(np.log10(grid.max() - grid.min()))

            vmax = np.around(2 * vmax / 10.**rounder) * 10.**rounder / 2.
            vmin = np.around(2 * vmin / 10.**rounder) * 10.**rounder / 2.

            if vmax == vmin:
                if vmax > 0:
                    vmin = 0
                else:
                    vmax = 0

        else:
            vmax = 1e-15
            vmin = -1e-15

    # Scale values if needed
    vmin *= scale
    vmax *= scale

    return vmin, vmax

def __adddots(ax: plt.axes,
              x: npt.ArrayLike,
              y: npt.ArrayLike,
              xlim: list,
              ylim: list):
    """
    Overlays circles onto an active 2D plot to show actual atomic positions.
    
    Parameters
    ----------
    ax : matplotlib.pyplot.axes
        The plotting axes to add the atomic positions to.
    x : array-like object
        List of x-coordinates for the atoms/dots.
    y : array-like object
        List of y-coordinates for the atoms/dots.
    xlim : list
        The [xmin, xmax] values associated with the plot.
    ylim : list
        The [ymin, ymax] values associated with the plot.
    """    
    # Set linewidth based on system dimensions
    syswidth = max([abs(xlim[1]-xlim[0]), abs(ylim[1]-ylim[0])])
    linewidth = 60. / syswidth

    # Add circles at each (x,y) position
    for xi, yi in zip(x, y):
        point = mpatches.Circle((xi, yi), 0.3, ls='solid', lw=linewidth)
        ax.add_artist(point)
        point.set_facecolor('none')
        point.set_edgecolor('k')
