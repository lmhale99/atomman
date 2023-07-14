# coding: utf-8
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://plotly.com/python/
try:
    import plotly.graph_objects as go
except ModuleNotFoundError:
    has_plotly = False
else:
    has_plotly = True

# atomman imports
from . import get_prop_values
from .. import System

__all__ = ['view_3d']

def view_3d(system: System,
            atomsize: float = 3,
            prop_name: str = 'atype',
            prop_index: Union[int, tuple, None] = None,
            prop_magnitude: bool = False,
            prop: Optional[npt.ArrayLike] = None,
            prop_unit: Optional[str] = None,
            colorscale: str = 'jet',
            width: int = 800,
            height: int = 800):
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
        The marker size to use to represent the atoms.  Note that in plotly,
        the absolute marker size is constant meaning that the relative size of
        the markers/atoms changes with box size and zooming.  Default value is 3.
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
    colorscale : str, optional
        The plotly colorscale to use for the atom colors.  Default value
        is 'jet'.
    width : int, optional
        The width in pixels to use for the plot size.  Default value is 800.
    height : int, optional
        The height in pixels to use for the plot size.  Default value is 800.
    """
    if has_plotly is False:
        raise ModuleNotFoundError('plotly is needed for this method')

    # Construct the cell box lines
    o = system.box.origin
    a = system.box.avect
    b = system.box.bvect
    c = system.box.cvect
    lines = np.array([o, o+a, o+a+c, o+a, o+a+b, o+a+b+c, o+a+b, o+b, o+b+c, o+b, o,
                      o+c, o+c+a, o+c+a+b, o+c+b, o+c])
    linesettings = {
        'color':'black',
    }
    figbox = go.Scatter3d(x=lines[:, 0], y=lines[:, 1], z=lines[:, 2], mode='lines',
                          line=linesettings)

    # Construct the atomic positions
    x = system.atoms.pos[:,0]
    y = system.atoms.pos[:,1]
    z = system.atoms.pos[:,2]
    p = get_prop_values(system, prop_name, prop_index, prop_magnitude, prop, prop_unit)[1]

    markersettings = {
        'size':atomsize,
        'color':p,
        'colorscale':colorscale,
    }
    figatoms = go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=markersettings)

    fig = go.Figure(data=[figatoms, figbox])

    fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False )
    fig.update_layout(
        width=width, height=height,
        showlegend=False,
        scene={
            'camera':{
                'up': {'x':0, 'y':0, 'z':1},
                'eye':{'x':0, 'y':2, 'z':0},
            },
            'aspectmode':'data'
        },
    )
    fig.show()
