import numpy as np
from scipy.interpolate import griddata

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import PatchCollection
          
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

def prettygrid(grid, xedges, yedges, cmap=cm.get_cmap('jet'), propname='', czero=True, scale=1):

    #This scales the color plot to be fancy
    #czero=True makes zero the halfway point of the scale
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
        
def a2cplot(system, value, index1=None, index2=None, plotbounds=None, bins=200, dots=True, czero=True, save=False, show=True):
    
    if save is False and show is False:
        print 'Figure not saved or displayed!'
    
    values, unit = system.atoms_prop(value)
    assert values is not None, 'atomic value not found!'
    
    name = str(value)
    if index1 is not None:
        if index2 is not None:
            name += '[' + str(index1) + '][' + str(index2) + ']'
            values = values[:, index1, index2]
        else:
            name += '[' + str(index1) + ']'
            values = values[:, index1]            
    
    assert (values.dtype == np.int or values.dtype == np.float) and values.shape == (system.natoms(), ), 'Value field is not a scalar number' 
    
    if plotbounds is None:
        plotbounds = np.array([[system.box('xlo')[0], system.box('xhi')[0]],
                               [system.box('ylo')[0], system.box('yhi')[0]],
                               [system.box('zlo')[0], system.box('zhi')[0]]])
    
    #Build arrays containing atomic information
    natoms = system.natoms()  
    pos = system.view['pos']
    x_check = np.all(pos[0] > plotbounds[0,0] - 5. and pos[0] < plotbounds[0,1] + 5.)
    y_check = np.all(pos[1] > plotbounds[1,0] - 5. and pos[1] < plotbounds[1,1] + 5.)
    z_check = np.all(pos[2] > plotbounds[2,0] and pos[2] < plotbounds[2,1])
    xyz_check = np.hstack((x_check[:, np.newaxis], y_check[:, np.newaxis], z_check[:, np.newaxis]))
    index_atoms = np.where(np.all(xyz_check), axis=1)[0]
    
    x = pos[index_atoms, 0]
    y = pos[index_atoms, 1]
    v = values[index_atoms]

    box = plotbounds[:2,:2]
    grid, xedges, yedges = grid_interpolate_2d(x, y, v, bins=bins, range=box)

    intsum = np.sum(grid)
    avsum = intsum / (bins-1) / (bins-1)
    intsum = intsum * (plotbounds[0,1]-plotbounds[0,0]) / (bins-1) * (plotbounds[1,1]-plotbounds[1,0]) / (bins-1)
        
    fig = prettygrid(grid, xedges, yedges, propname=name, czero=czero)
    if dots:
        adddots(x, y, xedges, yedges)
    if save:
        plt.savefig(name + '.png', dpi=800)
    if show==False:
        plt.close(fig)
    plt.show() 
 
    return intsum, avsum