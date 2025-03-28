# coding: utf-8
# Standard Python libraries
from re import A
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
from . import Strain
from ..plot import interpolate_contour

class DifferentialDisplacement():
    def __init__(self,
                 system0: System,
                 system1: System,
                 neighbors: Optional[NeighborList] = None,
                 cutoff: Optional[float] = None,
                 reference: int = 1):
        """
        Class initializer.  Calls solve if either neighbors or cutoff are given.
        
        Parameters
        ----------
        system0 : atomman.system
            The base/reference system to use.
        system1 : atomman.system
            The defect/current system to use.
        neighbors : atomman.NeighborList, optional
            The neighbor list to use.  
        cutoff : float, optional
            Cutoff distance for computing a neighbor list. If reference = 0, then system_0
            will be used to generate the list. If reference = 1, then system_1 will be
            used to generate the list.
        reference : int, optional
            Indicates which of the two systems should be used for the plotting
            reference: 0 or 1. If 0, then system_0's atomic positions will be
            used for the calculation and neighbors should be for system_0.  If
            1 (default), then system_1's atomic positions will be used
            for the calculation and neighbors should be for system_1.   
        """
        
        if neighbors is not None or cutoff is not None:
            self.solve(system0, system1, neighbors=neighbors, cutoff=cutoff, reference=reference)
        else:
            assert system0.natoms == system1.natoms
            
            self.__system0 = system0
            self.__system1 = system1
            self.reference = reference
            self.__neighbors = None
            self.__ddvectors = None
            self.__arrowcenters = None
            self.__arrowuvectors = None
    
    @property
    def reference(self) -> int:
        """int : Indicates which system (0 or 1) is used as the reference."""
        return self.__reference

    @reference.setter
    def reference(self, value: int):
        assert value == 0 or value == 1, 'reference must be 0 or 1'
        self.__reference = value    
    
    @property
    def system0(self) -> System:
        """atomman.System : The defect-free base system."""
        return self.__system0
    
    @property
    def system1(self) -> System:
        """atomman.System : The defect containing system."""
        return self.__system1
    
    @property
    def neighbors(self) -> NeighborList:
        """atomman.NeighborList : The list of neighbors identified for the reference system."""
        return self.__neighbors
    
    @property
    def ddvectors(self) -> Optional[np.ndarray]:
        """numpy.array or None : The computed differential displacement vectors."""
        return self.__ddvectors
    
    @property
    def arrowcenters(self) -> Optional[np.ndarray]:
        """numpy.array or None : The identified center positions for the ddvectors."""
        return self.__arrowcenters
    
    @property
    def arrowuvectors(self) -> Optional[np.ndarray]:
        """numpy.array or None : The unit vectors between all pairs of atoms for which the ddvectors have been computed."""
        return self.__arrowuvectors
    
    def solve(self,
              system0: Optional[System] = None,
              system1: Optional[System] = None,
              neighbors: Optional[NeighborList] = None,
              cutoff: Optional[float] = None,
              reference: Optional[int] = None):
        """
        Solves the differential displacement vectors.
        
        Parameters
        ----------
        system0 : atomman.system, optional
            The base/reference system to use.
        system1 : atomman.system, optional
            The defect/current system to use.
        neighbors : atomman.NeighborList, optional
            The neighbor list to use.  
        cutoff : float, optional
            Cutoff distance for computing a neighbor list. If reference = 0, then system_0
            will be used to generate the list. If reference = 1, then system_1 will be
            used to generate the list.
        reference : int, optional
            Indicates which of the two systems should be used for the plotting reference: 0 or 1.
            If 0, then system0's atomic positions will be used for the calculation and
            neighbors should be for system0.  If 1, then system1's atomic positions will be used
            for the calculation and neighbors should be for system1.  Default value is
            whatever was set when the object was initialized. 
        """
        # Handle parameters
        if system0 is not None:
            self.__system0 = system0
        else:
            system0 = self.system0
        if system1 is not None:
            self.__system1 = system1
        else:
            system1 = self.system1
        assert system0.natoms == system1.natoms
        
        if reference is None:
            reference = self.reference
        else:
            self.reference = reference
        if reference == 0:
            refsystem = system0
        else:
            refsystem = system1
        
        if neighbors is None:
            if cutoff is not None:
                self.__neighbors = neighbors = refsystem.neighborlist(cutoff=cutoff)
            else:
                if self.neighbors is not None:
                    neighbors = self.neighbors
                else:
                    raise ValueError('Either neighbors or cutoff must be given')
        else:
            self.__neighbors = neighbors
        
        all_ddvectors = []
        all_arrowcenters = []
        all_arrowuvectors = []
        
        # Loop over all atoms i in ref system
        for i in np.arange(refsystem.natoms):
            neighs = neighbors[i]
            if len(neighs) == 0:
                continue
            
            # Compute distance vectors between atom i and its neighbors for both systems
            dvectors0 = system0.dvect(int(i), neighs)
            dvectors1 = system1.dvect(int(i), neighs)
            if dvectors0.shape == (3,):
                dvectors0.shape = (1,3)
                dvectors1.shape = (1,3)

            # Compute differential displacement vectors
            ddvectors = dvectors1 - dvectors0

            # Compite center points and direction vectors
            if reference == 0:
                arrowcenters = system0.atoms.pos[i] + dvectors0 / 2
                arrowuvectors = dvectors0 / np.linalg.norm(dvectors0, axis=1)[:,np.newaxis]    
            else:    
                arrowcenters = system1.atoms.pos[i] + dvectors1 / 2
                arrowuvectors = dvectors1 / np.linalg.norm(dvectors1, axis=1)[:,np.newaxis]
                
            # Append calculation values to associated lists            
            all_ddvectors.append(ddvectors)
            all_arrowcenters.append(arrowcenters)
            all_arrowuvectors.append(arrowuvectors)
        
        # Save computed values to object properties 
        self.__ddvectors = np.concatenate(all_ddvectors)
        self.__arrowcenters = np.concatenate(all_arrowcenters)
        self.__arrowuvectors = np.concatenate(all_arrowuvectors)
        
    def plot(self,
             component: Union[str, npt.ArrayLike],
             ddmax: Optional[float],
             plotxaxis: Union[str, npt.ArrayLike] = 'x',
             plotyaxis: Union[str, npt.ArrayLike] = 'y',
             xlim: Optional[tuple] = None,
             ylim: Optional[tuple] = None,
             zlim: Optional[tuple] = None,
             arrowscale: float = 1,
             arrowwidth: float = 0.005,
             use0z: bool = False,
             atomcolor: Union[str, list, None] = None,
             atomcmap: Union[str, list, None] = None,
             atomsize: float = 0.5,
             figsize: int = 10,
             matplotlib_axes: Optional[plt.axes] = None
             ) -> Optional[plt.figure]:

        r"""
        Creates a matplotlib figure of a differential displacement map.  Atom
        positions are represented as circles, while the selected components of the
        differential displacement vectors are plotted as arrows.
        
        Parameters
        ----------
        component : str or array-like object
            Indicates the component(s) of the differential displacement to plot.
            Values of 'x', 'y', or 'z' will plot the component along that
            Cartesian direction.  A value of 'projection' will plot the
            differential displacement vectors as projected onto the plotting
            plane, thereby showing the two components perpendicular to the line
            direction.  If a 3D vector is given, then the component parallel to
            that direction will be used.
        ddmax : float or None
            The maximum differential displacement value allowed. Values will be
            kept between +-ddmax by wrapping values with larger absolute values
            around by adding/subtracting 2*ddmax. Typically, this is set to be
            \|b\|/2, but can be defect-specific. For instance, fcc a/2<110>
            dislocations and basal hcp dislocations are typically plotted with
            ddmax=\|b\|/4.  If set to None, then no wrapping is done.
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
        arrowscale : float, optional
            Scaling factor for the magnitude of the differential displacement
            arrows.  Default value is 1: no scaling, vectors are in units of length.
            For major components, this is often set such that the max differential
            displacement component after wrapping (see ddmax) is scaled to the
            distance between the atom pairs in the plot.  For minor components, this
            is often set to a large value simply to make the components visible.
        arrowwidth : float, optional
            Scaling factor to use for the width of the plotted arrows. Default value is
            0.005 = 1/200.
        use0z : bool, optional
            If False (default), the z coordinates from the reference system will be
            used for zlim and atomcmap colors. If True, the z coordinates will be
            used from system0 even if system1 is the reference system.
        atomcolor : str or list, optional
            Matplotlib color name(s) to use to display the atoms.  If str, that
            color will be assigned to all atypes.  If list, must give a color value
            or None for each atype.  Default value (None) will use cmap instead.
            Note: atomcolor and atomcmap can be used together as long as exactly
            one color or cmap is given for each unique atype.
        atomcmap : str or list, optional
            Matplotlib colormap name(s) to use to display the atoms.  Atoms will
            be colored based on their initial positions and scaled using zlim. If
            str, that cmap will be assigned to all atypes.  If list, must give a 
            cmap value or None for each atype.  Default value (None) will use 'hsv'
            cmap.  Note: atomcolor and atomcmap can be used together as long as
            exactly one color or cmap is given for each unique atype.
        atomsize : float, optional
            The circle radius size to use for the plotted atom positions in units of
            length.  Default value is 0.5.
        figsize : float or tuple, optional
            Specifies the size of the figure to create in inches.  If a single value
            is given, it will be used for the figure's width, and the height will be
            scaled based on the xlim and ylim values.  Alternatively, both the width
            and height can be set by passing a tuple of two values, but the plot will
            not be guaranteed to be "regular" with respect to length dimensions.
        matplotlib_axes : matplotlib.Axes.axes, optional
            An existing matplotlib axes object. If given, the differential displacement
            plot will be added to the specified axes of an existing figure.  This
            allows for subplots to be constructed.  Note that figsize will be ignored
            as the figure would have to be created beforehand and no automatic
            optimum scaling of the figure's dimensions will occur.
            
        Returns
        -------
        matplotlib.Figure
            The generated figure.  Not returned if matplotlib_axes is given.
        """

        ###################### Parameter and plot setup ########################

        # Interpret plot axis values
        plotxaxis = self.__plotaxisoptions(plotxaxis)
        plotyaxis = self.__plotaxisoptions(plotyaxis)

        # Build transformation matrix, T, from plot axes.
        T = axes_check([plotxaxis, plotyaxis, np.cross(plotxaxis, plotyaxis)])

        # Extract positions and transform using T        
        if self.reference == 0:        
            atompos = np.inner(self.system0.atoms.pos, T)
            refsystem = self.system0
        else:
            atompos = np.inner(self.system1.atoms.pos, T)
            refsystem = self.system1
            if use0z:
                pos0 = np.inner(self.system0.atoms.pos, T)
                atompos[:, 2] = pos0[:, 2]

        # Set default plot limits
        if xlim is None:
            xlim = (atompos[:, 0].min(), atompos[:, 0].max())
        if ylim is None:
            ylim = (atompos[:, 1].min(), atompos[:, 1].max())
        if zlim is None:
            zlim = (atompos[:, 2].min(), atompos[:, 2].max()) 

        # Define box for identifying only points inside
        plotbox = Box(xlo = xlim[0] - 5, xhi = xlim[1] + 5,
                      ylo = ylim[0] - 5, yhi = ylim[1] + 5,
                      zlo = zlim[0], zhi = zlim[1])
            
        # Set plot height if needed
        if isinstance(figsize, (int, float)):
            dx = xlim[1] - xlim[0]
            dy = ylim[1] - ylim[0]
            figsize = (figsize, figsize * dy / dx)

        # Initial plot setup and parameters
        if matplotlib_axes is None:
            fig = plt.figure(figsize=figsize, dpi=72)
            ax1 = fig.add_subplot(111)
        else:
            ax1 = matplotlib_axes
        ax1.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
        
        # Handle atomcolor and atomcmap values
        atomcolor, atomcmap = self.__atomcoloroptions(atomcolor, atomcmap)

        ######################## Add atom circles to plot ##############################
        
        # Loop over all atoms i in plotting box
        for i in np.arange(refsystem.natoms)[plotbox.inside(atompos)]:
            atype = refsystem.atoms.atype[i]
            atype_index = refsystem.atypes.index(atype)

            # Plot a circle for atom i
            if atomcmap[atype_index] is not None:
                color = atomcmap[atype_index]((atompos[i, 2] - zlim[0]) / (zlim[1] - zlim[0]))
            elif atomcolor[atype_index] is not None:
                color = atomcolor[atype_index]
            else:
                color = None
            if color is not None:
                ax1.add_patch(mpatches.Circle(atompos[i, :2], atomsize, fc=color, ec='k'))

        ######################## Arrow setup ##############################
        
        # Build arrows based on component
        arrowlengths, arrowcenters = self.__buildarrows(T, plotbox, component, ddmax)
        
        # Scale arrows
        arrowlengths = arrowscale * arrowlengths
        
        # Compute arrow widths based on lengths
        arrowwidths = arrowwidth * (arrowlengths[:,0]**2 + arrowlengths[:,1]**2)**0.5

        # Plot the arrows
        for center, length, width in zip(arrowcenters, arrowlengths, arrowwidths):
            if width > 1e-7:
                ax1.quiver(center[0], center[1], length[0], length[1],
                           pivot='middle', angles='xy', scale_units='xy',
                           scale=1, width=width, minshaft=2)
        
        if matplotlib_axes is None:
            return fig
    
    def __plotaxisoptions(self, plotaxis):
        """Internal method for handling plotxaxis and plotyaxis values"""
        
        # Give numeric values for str plot axis terms
        if plotaxis == 'x':
            plotaxis = [1.0, 0.0, 0.0]
        elif plotaxis == 'y':
            plotaxis = [0.0, 1.0, 0.0]
        elif plotaxis == 'z':
            plotaxis = [0.0, 0.0, 1.0]
        
        # Convert to numpy array
        return np.asarray(plotaxis, dtype=float)
    
    def __atomcoloroptions(self, atomcolor, atomcmap):
        """Internal method for handling atomcolor and atomcmap options"""
        
        # Identify number of atom types in the reference system
        if self.reference == 0:
            natypes = self.system0.natypes
        else:
            natypes = self.system1.natypes
        
        # Set default atomcmap
        if atomcolor is None and atomcmap is None:
            atomcmap = 'hsv'

        # Normalize atomcolor values
        if isinstance(atomcolor, str):
            atomcolor = [atomcolor for i in range(natypes)]
        elif atomcolor is None:
            atomcolor = [None for i in range(natypes)]
        else:
            atomcolor = list(atomcolor)
            if len(atomcolor) != natypes:
                raise ValueError('Invalid number of atomcolor values')

        # Normalize atomcmap values
        if isinstance(atomcmap, str):
            atomcmap = [atomcmap for i in range(natypes)]
        elif atomcmap is None:
            atomcmap = [None for i in range(natypes)]
        else:
            atomcmap = list(atomcmap)
            if len(atomcmap) != natypes:
                raise ValueError('Invalid number of atomcmap values')
        
        # Check that no atype has both atomcmap and atomcolor
        for color, cmap in zip(atomcolor, atomcmap):
            if color is not None and cmap is not None:
                raise ValueError('atomcmap and atomcolor cannot both be given for the same atype')
        
        # Convert atomcmap str values to color map objects
        for ic in range(natypes):
            if atomcmap[ic] is not None:    
                atomcmap[ic] = cm.get_cmap(atomcmap[ic])
        
        return atomcolor, atomcmap
    
    def __buildarrows(self, T, plotbox, component, ddmax):
        """Internal method for building the parameters for plotting the arrows"""
        
        # Manage component
        if isinstance(component, str):
            if component == 'x':
                component = np.array([1.0, 0.0, 0.0])
            elif component == 'y':
                component = np.array([0.0, 1.0, 0.0])
            elif component == 'z':
                component = np.array([0.0, 0.0, 1.0])
            elif component != 'projection':
                raise ValueError('Invalid component style: must be x, y, z, projection, or numpy array')
        else:
            component = np.asarray(component, dtype=float)
            assert component.shape == (3,), 'Invalid numeric component: must be a 3D vector'
            component = component / np.linalg.norm(component)
            
        # Transform arrow-related vectors
        ddvectors = np.inner(self.ddvectors, T)
        arrowcenters = np.inner(self.arrowcenters, T)
        arrowuvectors = np.inner(self.arrowuvectors, T)
                
        # Identify only vectors in ploting box
        inbounds = plotbox.inside(arrowcenters)
        ddvectors = ddvectors[inbounds]
        arrowcenters = arrowcenters[inbounds]
        arrowuvectors = arrowuvectors[inbounds]

        # Build arrows for the xy component option
        if isinstance(component, str) and component == 'projection':
            
            # Arrows are in-plane vector components
            ddcomponents = (ddvectors[:,0]**2 + ddvectors[:,1]**2)**0.5
            arrowuvectors = ddvectors[:, :2] / ddcomponents[:,np.newaxis]
            
            # Scheme for direction uniqueness (not sure what other projects use?)
            arrowuvectors[ddvectors[:, 2] > 0] *= -1
            
            # Normalize ddcomponents to be between +-ddmax
            if ddmax is not None and ddmax > 0:
                while True:
                    mask = ddcomponents > ddmax
                    if np.sum(mask) == 0:
                        break
                    ddcomponents[mask] -= 2 * ddmax
                while True:
                    mask = ddcomponents < -ddmax
                    if np.sum(mask) == 0:
                        break
                    ddcomponents[mask] += 2 * ddmax

            # Arrows have magnitude of ddcomponent
            arrowlengths = arrowuvectors * ddcomponents[:,np.newaxis]
        
        # Build arrows for vector components
        else:
            component = T.dot(component)
            ddcomponents = ddvectors.dot(component)

            # Normalize ddcomponents to be between +-ddmax
            if ddmax is not None and ddmax > 0:
                while True:
                    mask = ddcomponents > ddmax
                    if np.sum(mask) == 0:
                        break
                    ddcomponents[mask] -= 2 * ddmax
                while True:
                    mask = ddcomponents < -ddmax
                    if np.sum(mask) == 0:
                        break
                    ddcomponents[mask] += 2 * ddmax

            # Arrows have magnitude of ddcomponent and direction of uvectors
            arrowlengths = arrowuvectors * ddcomponents[:,np.newaxis]
        
        return arrowlengths, arrowcenters

    def plot_with_nye(self,
                      component: Union[str, npt.ArrayLike],
                      ddmax: Optional[float],
                      strain: Strain,
                      plotxaxis: Union[str, npt.ArrayLike] = 'x',
                      plotyaxis: Union[str, npt.ArrayLike] = 'y',
                      xlim: Optional[tuple] = None,
                      ylim: Optional[tuple] = None,
                      zlim: Optional[tuple] = None,
                      vlim: Optional[tuple] = None,
                      cmap: str = 'bwr',
                      arrowscale: float = 1,
                      arrowwidth: float = 0.005,
                      use0z: bool = False,
                      atomcolor: Union[str, list, None] = None,
                      atomcmap: Union[str, list, None] = None,
                      atomsize: float = 0.5,
                      figsize: int = 10,
                      xbins: int = 200,
                      ybins: int = 200,
                      colorbar: bool = True,
                      fill_value: float = np.nan,
                      matplotlib_axes: Optional[plt.axes] = None,
                      ) -> Tuple[float, plt.figure]:

        """
        Utility function for simple combined dd-Nye plots.  Note that this
        method is currently limited to components and plot axis values of
        'x', 'y', and 'z'.  This could be generalized in the future, if 
        there is interest...

        Parameters
        ----------
        component : str or array-like object
            Indicates the component(s) of the differential displacement to plot.
            Values of 'x', 'y', or 'z' will plot the component along that
            Cartesian direction.  A value of 'projection' will plot the
            differential displacement vectors as projected onto the plotting
            plane, thereby showing the two components perpendicular to the line
            direction.  If a 3D vector is given, then the component parallel to
            that direction will be used.
        ddmax : float or None
            The maximum differential displacement value allowed. Values will be
            kept between +-ddmax by wrapping values with larger absolute values
            around by adding/subtracting 2*ddmax. Typically, this is set to be
            |b|/2, but can be defect-specific. For instance, fcc a/2<110>
            dislocations and basal hcp dislocations are typically plotted with
            ddmax=|b|/4.  If set to None, then no wrapping is done.
        strain : atomman.defect.Strain
            A strain object computed for the system in question.  This will be used
            to compute the Nye tensor values that are included in the plot.
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
        vlim : tuple, optional
            Range limits for the Nye tensor contour plot.  If not given, will be
            determined by the range of the Nye tensor component values.
        cmap : str, optional
            The matplotlib colormap to use for the Nye tensor contour plot.
            Default value is 'bwr'.
        arrowscale : float, optional
            Scaling factor for the magnitude of the differential displacement
            arrows.  Default value is 1: no scaling, vectors are in units of length.
            For major components, this is often set such that the max differential
            displacement component after wrapping (see ddmax) is scaled to the
            distance between the atom pairs in the plot.  For minor components, this
            is often set to a large value simply to make the components visible.
        arrowwidth : float, optional
            Scaling factor to use for the width of the plotted arrows. Default value is
            0.005 = 1/200.
        use0z : bool, optional
            If False (default), the z coordinates from the reference system will be
            used for zlim and atomcmap colors. If True, the z coordinates will be
            used from system0 even if system1 is the reference system.
        atomcolor : str or list, optional
            Matplotlib color name(s) to use to display the atoms.  If str, that
            color will be assigned to all atypes.  If list, must give a color value
            or None for each atype.  Default value (None) will use cmap instead.
            Note: atomcolor and atomcmap can be used together as long as exactly
            one color or cmap is given for each unique atype.
        atomcmap : str or list, optional
            Matplotlib colormap name(s) to use to display the atoms.  Atoms will
            be colored based on their initial positions and scaled using zlim. If
            str, that cmap will be assigned to all atypes.  If list, must give a 
            cmap value or None for each atype.  Default value (None) will use 'hsv'
            cmap.  Note: atomcolor and atomcmap can be used together as long as
            exactly one color or cmap is given for each unique atype.
        atomsize : float, optional
            The circle radius size to use for the plotted atom positions in units of
            length.  Default value is 0.5.
        figsize : float or tuple, optional
            Specifies the size of the figure to create in inches.  If a single value
            is given, it will be used for the figure's width, and the height will be
            scaled based on the xlim and ylim values.  Alternatively, both the width
            and height can be set by passing a tuple of two values, but the plot will
            not be guaranteed to be "regular" with respect to length dimensions.
        xbins : int, optional
            Specifies the number of interpolation bins to use along the plotting 
            x-axis.  Default value is 200.
        ybins : int, optional
            Specifies the number of interpolation bins to use along the plotting
            y-axis.  Default value is 200.
        colorbar: bool, optional
            If True (default) a colorbar will be added to the plot.
        fill_value: float, optional
            Value used to fill in for grid points failed to interpolate in the fit.
            If not given, then the default is np.nan, which may cause an error in
            plotting for too narrow xlim and ylim settings.
        matplotlib_axes : matplotlib.Axes.axes, optional
            An existing matplotlib axes object. If given, the differential displacement
            plot will be added to the specified axes of an existing figure.  This
            allows for subplots to be constructed.  Note that figsize will be ignored
            as the figure would have to be created beforehand and no automatic
            optimum scaling of the figure's dimensions will occur.

        Returns
        -------
        float
            The integer sum of the Nye tensor over the plotted area.
        matplotlib.Figure
            The generated figure.

        """
        # Identify Nye tensor components based on component and plot axes
        try: 
            component_index = ['x', 'y', 'z'].index(component)
        except:
            raise ValueError('component is currently limited to values of "x", "y" or "z"') 
        try:
            x_index = ['x', 'y', 'z'].index(plotxaxis)
        except:
            raise ValueError('plotxaxis is currently limited to values of "x", "y" or "z"')
        try:
            y_index = ['x', 'y', 'z'].index(plotyaxis)
        except:
            raise ValueError('plotyaxis is currently limited to values of "x", "y" or "z"')
        line_index = 3 - (x_index + y_index)
        prop_index = [line_index, component_index] 
        
        # Generate dd plot
        if matplotlib_axes is None:
            fig = self.plot(component, ddmax, plotxaxis=plotxaxis, plotyaxis=plotyaxis,
                            xlim=xlim, ylim=ylim, zlim=zlim, arrowscale=arrowscale,
                            arrowwidth=arrowwidth, use0z=use0z, atomcolor=atomcolor,
                            atomcmap=atomcmap, atomsize=atomsize, figsize=figsize)
            matplotlib_axes = fig.axes[0]
        else:
            self.plot(component, ddmax, plotxaxis=plotxaxis, plotyaxis=plotyaxis,
                      xlim=xlim, ylim=ylim, zlim=zlim, arrowscale=arrowscale,
                      arrowwidth=arrowwidth, use0z=use0z, atomcolor=atomcolor,
                      atomcmap=atomcmap, atomsize=atomsize, figsize=figsize,
                      matplotlib_axes=matplotlib_axes)
            fig = None

        # Add Nye tensor surface plot
        intsum = interpolate_contour(self.system1, 'nye', prop_index=prop_index, prop=strain.nye,
                            plotxaxis=plotxaxis, plotyaxis=plotyaxis,
                            xlim=xlim, ylim=ylim, zlim=zlim, vlim=vlim, cmap=cmap,
                            xbins=xbins, ybins=ybins, fill_value=fill_value, 
                            matplotlib_axes=matplotlib_axes,
                            dots=False, colorbar=colorbar, title=False)[0]
        
        if fig is None:
            return intsum
        else:
            return intsum, fig