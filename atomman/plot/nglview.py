# coding: utf-8
from typing import Optional, Union

# http://www.numpy.org/
import numpy.typing as npt

# Switch for if nglview version has show_atomman or not
try:
    from nglview.show import show_atomman
except (ModuleNotFoundError, ImportError):
    from .nglview_classes import show_atomman, has_nglview
else:
    has_nglview = True

# atomman imports
from . import get_prop_values, values_to_hexcolors
from .. import System

def view_3d(system: System,
            atomsize: float = 1,
            prop_name: str = 'atype',
            prop_index: Union[int, tuple, None] = None,
            prop_magnitude: bool = False,
            prop: Optional[npt.ArrayLike] = None,
            prop_unit: Optional[str] = None,
            cmap: str = 'jet',
            cmin: Optional[float] = None,
            cmax: Optional[float] = None,
            #width: int = 800,
            #height: int = 800,
            ):
    """
    Visualize an atomic configuration in 3D using plotly.  Useful for quick
    visualizations during Jupyter Notebook sessions.

    Note: This is a new feature and supported parameters, parameter names, and
    even the name of the function may be subject to change.

    Parameters
    ----------
    system : atomman.System
        The atomic system to generate the visualization for.
    atomsize : float
        The radius size for the atomic spheres.  Default value is 1.
    prop_name : str
        The name of the per-atom property stored in system that you want to
        color the atoms by.  Default value is 'atype'.
    prop_index : int or tuple, optional
        If prop_name is a multi-dimensional array, prop_index allows for a
        specific index to be selected to plot.  Cannot be used with
        prop_magnitude as they are competing options.
    prop_magnitude : bool, optional
        If prop_name is a per-atom vector, setting prop_magnitude as True will
        use the magnitude of the prop_name vector.  Cannot be used with
        prop_index as they are competing options.
    prop : array-like object, optional
        Allows for the property values to color the atoms by to be given
        directly rather than stored as a per-atom property of system.
    prop_unit : str or None, optional
        The units to use for the property value being plotted.  Default value
        is None, in which no unit conversion is applied.
    cmap : str, optional
        The plotly colorscale to use for the atom colors.  Default value
        is 'jet'.
    cmin : float or None, optional
        The minimum limit value to use for the color scale.  If None (default)
        then the minimum value in values will be used.
    cmax : float or None, optional
        The maximum limit value to use for the color scale.  If None (default)
        then the maximum value in values will be used.
    width : int, optional
        The width in pixels to use for the plot size.  Default value is 800.
    height : int, optional
        The height in pixels to use for the plot size.  Default value is 800.
    """
    if has_nglview is False:
        raise ModuleNotFoundError('nglview is needed for this method')

    # Get property values based on inputs
    values = get_prop_values(system, prop_name, prop_index, prop_magnitude, prop, prop_unit)[1]

    # Convert input values to hex colors
    colors = values_to_hexcolors(values, cmap, cmin, cmax)

    view = show_atomman(system)

    #view.clear_representations()
    #for i, color in enumerate(colors):
    #    view.add_representation("ball+stick", selection=f"[{i+1}]", color=color)
    view.add_unitcell()

    return view
