'''
Find paths defining 3km across- and along-shelf boxes for examining 
coastal connectivity.
'''

import numpy as np
import tracpy
import tracpy.plotting
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt

grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)

# To click to get the points:
# tracpy.plotting.background(grid)
# plt.plot(grid['xpsi'], grid['ypsi'], 'k', grid['xpsi'].T, grid['ypsi'].T, 'k')
# ind = grid['mask'].astype(bool)
# plt.plot(grid['xr'][ind], grid['yr'][ind], 'bs')
# pts = ginput(timeout=0, n=0)
# np.savez('calcs/coastpts.npz', lon=lon, lat=lat)

# Read in previously-clicked rho points defining the inner edge of 
# the inner shelf (excludes the bays)
d = np.load('calcs/coastpts.npz')
lon = d['lon']; lat = d['lat']
d.close()
x, y = grid['basemap'](lon, lat)

## Add in more points
# Calculate distance along-shore for each line
dd = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
d = np.zeros(x.size)
d[1:] = np.cumsum(dd)
while dd.max()>100:
    xnew = np.empty(x.size*2-1)
    xnew[1::2] = .5*(x[1:]+x[:-1])
    xnew[::2] = x
    ynew = np.empty(y.size*2-1)
    ynew[1::2] = .5*(y[1:]+y[:-1])
    ynew[::2] = y
    x = xnew.copy(); y = ynew.copy()
    dd = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
    d = np.zeros(x.size)
    d[1:] = np.cumsum(dd)

# spacing at 3 km along the coast
d3 = np.arange(d.min(), d.max(), 3000)
x3 = np.interp(d3, d, x)
y3 = np.interp(d3, d, y)

# Calculate normal vector
dx = x3[1:] - x3[:-1]
dy = y3[1:] - y3[:-1]
mag = np.sqrt(dx**2 + dy**2)
xn = x3[:-1] + (dy/mag)*3000
yn = y3[:-1] - (dx/mag)*3000

# put the 4 parts together to make path for each box along-shore
paths = []
for i in xrange(xn.size-1):
    verts = [(x3[i+1],y3[i+1]), (xn[i+1],yn[i+1]), 
            (xn[i],yn[i]), (x3[i],y3[i])]
    paths.append(Path(verts))

# Plot up
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid, ax=ax)
plt.plot(grid['xpsi'], grid['ypsi'], 'k', grid['xpsi'].T, grid['ypsi'].T, 'k')
ind = grid['mask'].astype(bool)
plt.plot(grid['xr'][ind], grid['yr'][ind], 'bs')
for path in paths:
    patch = patches.PathPatch(path, facecolor='orange', lw=2, zorder=10)
    ax.add_patch(patch)