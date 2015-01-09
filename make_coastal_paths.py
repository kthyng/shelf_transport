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

ds = 5 # spacing at 3 km along the coast
d3 = np.arange(d.min(), d.max(), 1000*ds)
x3 = np.interp(d3, d, x)
y3 = np.interp(d3, d, y)

# Calculate normal vector
dx = x3[1:] - x3[:-1]
dy = y3[1:] - y3[:-1]
mag = np.sqrt(dx**2 + dy**2)
xn = x3[:-1] + (dy/mag)*1000*ds
yn = y3[:-1] - (dx/mag)*1000*ds

lonn, latn = grid['basemap'](xn, yn, inverse=True)
lon3, lat3 = grid['basemap'](x3, y3, inverse=True)
xg3, yg3, _ = tracpy.tools.interpolate2d(x3, y3, grid, 'd_xy2ij')
xgn, ygn, _ = tracpy.tools.interpolate2d(xn, yn, grid, 'd_xy2ij')

# put the 4 parts together to make path for each box along-shore
paths = []; pathsxy = []; pathsg = []
for i in xrange(xn.size-1):
    vertsg = [(xg3[i+1],yg3[i+1]), (xgn[i+1],ygn[i+1]), 
            (xgn[i],ygn[i]), (xg3[i],yg3[i])]
    vertsxy = [(x3[i+1],y3[i+1]), (xn[i+1],yn[i+1]), 
            (xn[i],yn[i]), (x3[i],y3[i])]
    verts = [(lon3[i+1],lat3[i+1]), (lonn[i+1],latn[i+1]), 
            (lonn[i],latn[i]), (lon3[i],lat3[i])]
    pathsg.append(Path(vertsg))
    pathsxy.append(Path(vertsxy))
    paths.append(Path(verts))

# make outer path of all boxes
verts_outer = np.vstack((np.vstack((xg3,yg3)).T, 
                        np.vstack((xgn[::-1],ygn[::-1])).T,
                        [xg3[0],yg3[0]]))
outerpath = Path(verts_outer)

np.savez('calcs/coastpaths.npz', paths=paths, pathsg=pathsg, 
            pathsxy=pathsxy, outerpathg=outerpath)

# Plot up
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid, ax=ax)
plt.plot(grid['xpsi'], grid['ypsi'], 'k', grid['xpsi'].T, grid['ypsi'].T, 'k')
ind = grid['mask'].astype(bool)
plt.plot(grid['xr'][ind], grid['yr'][ind], 'bs')
for path in pathsxy:
    patch = patches.PathPatch(path, facecolor='orange', lw=2, zorder=10)
    ax.add_patch(patch)