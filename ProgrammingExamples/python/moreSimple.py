# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 12:53:32 2012

@author: alfaceor
"""

import numpy as np
import matplotlib.pyplot as plt
x,y,z = np.loadtxt("datatest.dat",usecols=(5,4,3)).T
xbinsize=0.5
ybinsize=0.1
xi = np.arange(x.min(),x.max()+xbinsize,xbinsize)
yi = np.arange(y.min(),y.max()+ybinsize,ybinsize)

xim,yim = np.meshgrid(xi,yi)

# colorbar stuff
# palette = plt.matplotlib.colors.Colormap()
#palette = plt.matplotlib.colors.LinearSegmentedColormap('jet3',plt.cm.datad['jet'],2048)
#palette.set_under(alpha=0.0)

grid = np.zeros(xim.shape,dtype=xim.dtype)
bins = np.copy(grid)

nrow, ncol = grid.shape

# Not a number array for the white color in the plot
for row in range(nrow):
    for col in range(ncol):
        # TODO: se tiene que testar si existe un
        grid[row, col]=np.nan

        xc = xim[row, col]    # x coordinate.
        yc = yim[row, col]    # y coordinate.

        # find the position that xc and yc correspond to.
        posx = np.abs(x - xc)
        posy = np.abs(y - yc)
        ibin = np.logical_and(posx < xbinsize/2., posy < ybinsize/2.)
        ind  = np.where(ibin == True)[0]
        # Allocate the values in corresponding grids
        zsum=0
        if len(ind) != 0:
            for val in ind:
                zsum=zsum+z[val]
            grid[row,col] = zsum/len(ind)
        else:
            grid[row,col] = np.nan

# ahora tengo q colocar los valores adecuados
extent = (x.min(), x.max(), y.min(), y.max()) # extent of the plot
#plt.contourf(xim,yim,grid)
plt.imshow(grid, extent=extent, cmap=None, origin='lower', vmin=z.min(), vmax=z.max(), aspect='auto', interpolation='bilinear')
plt.colorbar()
plt.scatter(x,y,marker='o',c='b',s=1)
plt.show()
#plt.scatter(x,y)
