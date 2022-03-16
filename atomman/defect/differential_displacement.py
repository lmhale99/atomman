# coding: utf-8
# Standard Python libraries
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://matplotlib.org/
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

# atomman imports
from ..tools import axes_check
from .. import Box, System, NeighborList

def differential_displacement(system_0: System,
                              system_1: System,
                              burgers: npt.ArrayLike,
                              plotxaxis: Union[str, npt.ArrayLike] = 'x',
                              plotyaxis: Union[str, npt.ArrayLike] = 'y',
                              xlim: Optional[tuple] = None,
                              ylim: Optional[tuple] = None,
                              zlim: Optional[tuple] = None,
                              neighbors: Optional[NeighborList] = None,
                              cutoff: Optional[float] = None,
                              component: str = 'standard',
                              axes: Optional[npt.ArrayLike] = None,
                              plot_scale: float = 1,
                              atom_color: Union[str, list, None] = None,
                              atom_cmap: Union[str, list, None] = None,
                              display_final_pos: bool = False,
                              return_data: bool = False,
                              matplotlib_axes: Optional[plt.axes] = None
                              ) -> Union[plt.figure, dict,
                                         Tuple[plt.figure, dict], None]:
    """
    Generates a differential displacement plot for characterizing dislocations.
    
    Parameters
    ----------
    system_0 : atomman.system
        The base/reference system to use.
    system_1 : atomman.system
        The defect/current system to use.
    burgers : array-like object
        The dislocation's Burgers vector.
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
        in the specified length_unit.  The optimum zlim should encompass only
        a single periodic slice.  If not given, then the limits are set
        based on min and max atomic coordinates along the axis.
    neighbors : atomman.NeighborList, optional
        The neighbor list associated with system_0 to use.  Either neighbors
        or cutoff must be given, or system_0 must have a neighbors attribute.
    cutoff : float, optional
        Cutoff distance for computing a neighbor list for system_0.  Either
        neighbors or cutoff must be given, or system_0 have a neighbors
        attribute.
    component : str, optional
        Indicates the style of the calculation to use.  'standard' (default)
        plots the differential displacements between atoms in the Burgers
        direction.  'xy' plots the differential displacements within the xy
        plotting plane.  This is useful for screw dislocations with localized
        non-screw components.
    axes : array-like object, optional
        3x3 transformation array to apply to the Burgers vector to make it
        correspond to the system's orientation.
    plot_scale : float, optional
        Scalar for multiplying the magnitude of the differential displacement
        arrows.  Default value is 1 (no scaling).
    atom_color : str or list, optional
        Matplotlib color name(s) to use to display the atoms.  If str, that
        color will be assigned to all atypes.  If list, must give a color value
        or None for each atype.  Default value (None) will use cmap instead.
        Note: atom_color and atom_cmap can be used together as long as exactly
        one color or cmap is given for each unique atype.
    atom_cmap : str or list, optional
        Matplotlib colormap name(s) to use to display the atoms.  Atoms will
        be colored based on their initial positions and scaled using zlim. If
        str, that cmap will be assigned to all atypes.  If list, must give a 
        cmap value or None for each atype.  Default value (None) will use 'hsv'
        cmap.  Note: atom_color and atom_cmap can be used together as long as
        exactly one color or cmap is given for each unique atype.
    display_final_pos : bool, optional
        Flag to display positions of atoms and arrows relative to final
        configuration (system_1) rather than initial configuration (system_0).
        Note that this does not affect the atom's cmap color as the initial
        plotzaxis is always used.
    return_data : bool, optional
        If True, will return a dict containing the differential displacement
        vectors and vector positions.  Default is False.  Note: returned values
        are oriented relative to the plotting axes.
    matplotlib_axes : matplotlib.Axes.axes, optional
        An existing matplotlib axes object. If given, the differential displacement
        plot will be added to the specified axes of an existing figure.  This
        allows for subplots to be constructed.

    Returns
    -------
    fig : matplotlib.figure
        The generated figure. Not returned if matplotlib_axes is given.
    data : dict
        Contains differential displacement vectors and arrow plotting information.
        Returned if return_data is True.
    """
    # Transform burgers using axes
    if axes is not None:
        T = axes_check(axes)
        burgers = T.dot(burgers)
    
    # Neighbor list setup
    if neighbors is not None:
        assert cutoff is None, 'neighbors and cutoff cannot both be given'
    elif cutoff is not None:
        neighbors = NeighborList(system=system_0, cutoff=cutoff)
    elif hasattr(system_0, 'neighbors'):
        neighbors = system_0.neighbors
    else:
        raise ValueError('neighbors or cutoff is required')
    
    # Give numeric values for str plot axis terms
    if plotxaxis == 'x':
        plotxaxis = [1.0,0.0,0.0]
    elif plotxaxis == 'y':
        plotxaxis = [0.0,1.0,0.0]
    elif plotxaxis == 'z':
        plotxaxis = [0.0,0.0,1.0]
    
    if plotyaxis == 'x':
        plotyaxis = [1.0,0.0,0.0]
    elif plotyaxis == 'y':
        plotyaxis = [0.0,1.0,0.0]
    elif plotyaxis == 'z':
        plotyaxis = [0.0,0.0,1.0]
    
    # Build transformation matrix, T, from plot axes.
    plotxaxis = np.asarray(plotxaxis, dtype='float64')
    plotyaxis = np.asarray(plotyaxis, dtype='float64')
    T = axes_check([plotxaxis, plotyaxis, np.cross(plotxaxis, plotyaxis)])
    
    # Extract positions and transform using T
    pos_0 = np.inner(system_0.atoms.pos, T)
    pos_1 = np.inner(system_1.atoms.pos, T)
    
    #  Select initial or final reference configurations
    if display_final_pos:
        pos_ref = pos_1
        system_ref = system_1
    else:
        pos_ref = pos_0
        system_ref = system_0

    # Transform burgers using plot axes and separate magnitude and direction
    burgers = T.dot(burgers)
    burgers_mag = np.linalg.norm(burgers)
    burgers_uvect = burgers / burgers_mag
    
    # Set plot limits
    if xlim is None:
        xlim = (pos_ref[:, 0].min(), pos_ref[:, 0].max())
    if ylim is None:
        ylim = (pos_ref[:, 1].min(), pos_ref[:, 1].max())
    if zlim is None:
        zlim = (pos_0[:, 2].min(), pos_0[:, 2].max()) # pylint: disable=invalid-sequence-index
    
    # Identify only atoms in xlim, ylim, zlim
    in_bounds = ((pos_ref[:,0] > xlim[0] - 5.) &
                 (pos_ref[:,0] < xlim[1] + 5.) &
                 (pos_ref[:,1] > ylim[0] - 5.) &
                 (pos_ref[:,1] < ylim[1] + 5.) &
                 (pos_0[:,2] > zlim[0]) & # pylint: disable=invalid-sequence-index
                 (pos_0[:,2] < zlim[1])) # pylint: disable=invalid-sequence-index

    # if matplotlib_axes are passed do not create new
    if matplotlib_axes is not None:
        ax1 = matplotlib_axes
    else:
        # Initial plot setup and parameters
        fig, ax1, = plt.subplots(1, 1, squeeze=True, figsize=(7,7), dpi=72)

    ax1.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
    atom_circle_radius = burgers_mag / 10
    arrow_width_scale = 1. / 200.
    
    # Set default atom_cmap
    if atom_color is None and atom_cmap is None:
        atom_cmap = 'hsv'
    
    # Transform single color/cmap values to lists
    if isinstance(atom_cmap, str):
        if atom_color is not None:
            raise TypeError('atom_cmap and atom_color cannot be str if both are given')
        atom_cmap = [atom_cmap for i in range(system_ref.natypes)]
        atom_color = [None for i in range(system_ref.natypes)]
    elif isinstance(atom_color, str):
        if atom_cmap is not None:
            raise TypeError('atom_cmap and atom_color cannot be str if both are given')
        atom_color = [atom_color for i in range(system_ref.natypes)]
        atom_cmap = [None for i in range(system_ref.natypes)]
    else:
        atom_color = list(atom_color)
        atom_cmap = list(atom_cmap)
    
    # Check atom_color, atom_cmap list compatibility
    if len(atom_cmap) != system_ref.natypes:
        raise ValueError('Invalid number of atom_cmap values')
    if len(atom_color) != system_ref.natypes:
        raise ValueError('Invalid number of atom_color values')
    for ic in range(len(atom_cmap)):
        if atom_cmap[ic] is not None:
            if atom_color[ic] is not None:
                raise ValueError('atom_cmap and atom_color cannot both be given for same atype')
            atom_cmap[ic] = cm.get_cmap(atom_cmap[ic])

    data = {}
    data['centers'] = []
    data['vectors'] = []

    # Loop over all atoms i in plot range
    for i in np.arange(system_ref.natoms)[in_bounds]:
        
        atype = system_ref.atoms.atype[i]
        atype_index = system_ref.atypes.index(atype)
        
        # Plot a circle for atom i
        if atom_cmap[atype_index] is not None:
            color = atom_cmap[atype_index]((pos_0[i, 2] - zlim[0]) / (zlim[1] - zlim[0])) # pylint: disable=invalid-sequence-index
        elif atom_color[atype_index] is not None:
            color = atom_color[atype_index]
        else:
            color = None
        if color is not None:
            ax1.add_patch(mpatches.Circle(pos_ref[i, :2], atom_circle_radius, fc=color, ec='k'))

        # Compute distance vectors between atom i and its neighbors for both systems
        dvectors_0 = np.inner(system_0.dvect(int(i), neighbors[i]), T)
        dvectors_1 = np.inner(system_1.dvect(int(i), neighbors[i]), T)
    
        # Compute differential displacement vectors
        dd_vectors = dvectors_1 - dvectors_0
        
        if display_final_pos:
            # Compute center point positions for the vectors
            arrow_centers = pos_1[i] + dvectors_1 / 2
        else :
            arrow_centers = pos_0[i] + dvectors_0 / 2
        
        if return_data:
            data['vectors'].append(dd_vectors)
            data['centers'].append(arrow_centers)

        # Plot standard differential displacement component
        if component == 'standard':
            # Compute unit distance vectors
            if display_final_pos:
                uvectors = dvectors_1 / np.linalg.norm(dvectors_1, axis=1)[:,np.newaxis]
            else:
                uvectors = dvectors_0 / np.linalg.norm(dvectors_0, axis=1)[:,np.newaxis]
            
            # Compute component of the dd_vector parallel to the burgers vector
            dd_components = dd_vectors.dot(burgers_uvect)
            dd_components[dd_components > burgers_mag / 2] -= burgers_mag
            dd_components[dd_components < -burgers_mag / 2] += burgers_mag
            
            # Scale arrow lengths and vectors
            arrow_lengths = uvectors * dd_components[:,np.newaxis] * plot_scale
            arrow_widths = arrow_width_scale * dd_components * plot_scale
            
            # Plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths,
                                             arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1],
                               pivot='middle', angles='xy', scale_units='xy',
                               scale=1, width=width)

        
        # Plot xy differential displacement component
        elif component == 'xy':
            # Scale arrow lengths and vectors
            arrow_lengths = dd_vectors[:, :2] * plot_scale
            arrow_lengths[dd_vectors[:, 2] > 0] *= -1
            arrow_widths = arrow_width_scale * (arrow_lengths[:,0]**2
                                              + arrow_lengths[:,1]**2)**0.5
            
            # Plot the arrows
            for center, length, width in zip(arrow_centers, arrow_lengths, arrow_widths):
                if width > 1e-7:
                    ax1.quiver(center[0], center[1], length[0], length[1], width=width,
                               pivot='middle', angles='xy', scale_units='xy', scale=1)
    
    returns = []
    if matplotlib_axes is None:
        returns.append(fig)
    
    if return_data:
        data['centers'] = np.concatenate(data['centers'])
        data['vectors'] = np.concatenate(data['vectors'])
        returns.append(data)

    if len(returns) == 1:
        return returns[0]
    
    elif len(returns) > 1:
        return tuple(returns)