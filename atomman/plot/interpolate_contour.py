import numpy as np
from scipy.interpolate import griddata

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import PatchCollection

import atomman.unitconvert as uc

def interpolate_contour(system, property, index=None, plotbounds=None, bins=200, dots=True, czero=True, save=False, show=True, length_unit='angstrom', property_unit=None, cmap='jet'):
    """Creates a contour plot of a system's per-atom properties by interpolating between all of the atoms."""
    
    if save is False and show is False:
        print 'Figure not saved or displayed!'
    
    values = system.atoms_prop(key=property)
    assert values is not None, 'atomic property not found!'
    
    name = str(property)
    if index is not None:
        if isinstance(index, (int, long)):
            index = [index]
        else:
            index = list(index)
        for i in index:
            name += '[' + str(i+1) + ']'
        values = values[[Ellipsis] + index]
    
    assert values.shape == (system.natoms, ), 'Identified property[index] must be a per-atom scalar' 
    
    if plotbounds is None:
        plotbounds = np.array([[system.box.xlo, system.box.xhi],
                               [system.box.ylo, system.box.yhi],
                               [system.box.zlo, system.box.zhi]])
        plotbounds = uc.get_in_units(plotbounds, length_unit)
        
    #Build arrays containing atomic information
    pos = uc.get_in_units(system.atoms_prop(key='pos'), length_unit)
    
    #identify only atoms in plotbounds
    in_bounds = ((pos[:,0] > plotbounds[0,0] - 5.) &
                 (pos[:,0] < plotbounds[0,1] + 5.) &
                 (pos[:,1] > plotbounds[1,0] - 5.) &
                 (pos[:,1] < plotbounds[1,1] + 5.) &
                 (pos[:,2] > plotbounds[2,0]) &
                 (pos[:,2] < plotbounds[2,1]))
    
    x = pos[in_bounds, 0]
    y = pos[in_bounds, 1]
    v = values[in_bounds]

    box = plotbounds[:2,:2]
    grid, xedges, yedges = grid_interpolate_2d(x, y, v, bins=bins, range=box)

    intsum = np.sum(grid)
    avsum = intsum / (bins-1) / (bins-1)
    intsum = intsum * (plotbounds[0,1]-plotbounds[0,0]) / (bins-1) * (plotbounds[1,1]-plotbounds[1,0]) / (bins-1)
        
    fig = prettygrid(grid, xedges, yedges, cmap=cmap, propname=name, czero=czero)
    if dots:
        adddots(x, y, xedges, yedges)
    if save:
        plt.savefig(name + '.png', dpi=800)
    if show==False:
        plt.close(fig)
    plt.show() 
 
    return intsum, avsum
    
def grid_interpolate_2d(x, y, v, bins=50, range=None):
    #Handle range and bins options
    if range is None:
        range=[[x.min() ,x.max()], [y.min(), y.max()]]
    if isinstance(bins, int):
        xi = np.linspace(range[0][0], range[0][1], num=bins)
        yi = np.linspace(range[1][0], range[1][1], num=bins)
    elif isinstance(bins, list) or isinstance(bob, tuple) or isinstance(bins, np.ndarray):
        if len(bins)==2 and isinstance(bins[0], int) and isinstance(bins[1], int):
            xi = np.linspace(range[0][0], range[0][1], num=bins[0])
            yi = np.linspace(range[1][0], range[1][1], num=bins[1])
    
    #Create continuum grid from atomic information
    x0, y0 = np.meshgrid(xi, yi)
    
    grid = griddata((x, y), v, (x0, y0))
    
    return grid, range[0], range[1]

def prettygrid(grid, xedges, yedges, cmap='jet', propname='', czero=True, scale=1):

    #This scales the color plot to be fancy
    #czero=True makes zero the halfway point of the scale
    
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)
    
    if czero:
        vmax = abs(grid).max()
        if vmax != 0.0:
            vrounder = np.floor(np.log10(vmax))
            vmax = np.around(2 * vmax / 10.**vrounder) * 10.**vrounder / 2.
        else:
            vmax = 1e-15
        vmin = -vmax
    else:
        vmax = grid.max()
        vmin = grid.min()
        if abs(grid).max() != 0.0:
            vrounder = np.floor(np.log10(abs(grid).max()))
            vmax = np.around(2 * vmax / 10.**vrounder) * 10.**vrounder / 2.
            vmin = np.around(2 * vmin / 10.**vrounder) * 10.**vrounder / 2.
            if vmax == vmin:
                if vmax > 0:
                    vmin = 0
                else:
                    vmax = 0
        else:
            vmax = 1e-15
            vmin = -1e-15
    vmin*=scale
    vmax*=scale
    
    vticks = np.linspace(vmin, vmax, 11, endpoint=True)

    #Plot figure on screen
    fig = plt.figure(figsize=(7.7, 7), dpi=72)
    plt.imshow(grid, extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]], 
               cmap=cmap, norm=plt.Normalize(vmax=vmax, vmin=vmin))
    
    #Make colorbar pretty
    cbar = plt.colorbar(fraction=0.05, pad = 0.05, ticks=vticks)
    cbar.ax.tick_params(labelsize=15)
   
    #Make axis values and title pretty
    plt.xlim(xedges[0], xedges[-1])
    plt.xticks(size=15)
    plt.ylim(yedges[0], yedges[-1])
    plt.yticks(size=15)
    plt.title(propname, size=30)
    return fig
    
def adddots(x, y, xedges, yedges):
    points = []
    ax = plt.gca()
    syswidth = max([abs(xedges[-1]-xedges[0]), abs(yedges[-1]-yedges[0])])
    linewidth = 60. / syswidth
    for i in xrange(len(x)):
        point = mpatches.Circle((x[i], y[i]), 0.3, ls='solid', lw=linewidth)
        ax.add_artist(point)
        point.set_facecolor('none')
        point.set_edgecolor('k')