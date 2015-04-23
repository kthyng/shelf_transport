'''
Plot some sample drifters that reach the coastline in winter and summer.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
import tracpy.calcs
from datetime import datetime, timedelta
import glob
import op
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar
import matplotlib.patches as Patches
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from collections import defaultdict
from matplotlib.path import Path
import matplotlib.patches as patches

mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


grid = tracpy.inout.readgrid('grid.nc', usebasemap=True)

verts_outerxy = np.load('calcs/coastpaths.npz')['outerpathxy'].item()
pathsxy = np.load('calcs/coastpaths.npz')['pathsxy']


# Plot drifters from a winter simulation
dd = 5000 # decimate the drifters
dcalcsW = np.load('calcs/coastconn/likelihood/hist-2008-01.npz')
dtracksW = netCDF.Dataset('tracks/2008-01-01T00gc.nc')
xg = dtracksW.variables['xg'][:]; yg = dtracksW.variables['yg'][:]
# xg = xg[::dd,:]; yg = yg[::dd,:]
ind = xg==-1

xpW, ypW, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xpW[ind] = np.nan; ypW[ind] = np.nan

# drifter ids
ids = dcalcsW['ids'].item()
idsuse = []
[idsuse.extend(ids[0,5,i]) for i in xrange(342)]
idsuseW = set(idsuse) # unique only
# set up a set of indices with respect to all drifters with 
# coastal ones as True
indsW = np.zeros(xpW.shape[0]).astype(bool)
indsW[list(idsuseW)] = True



# Plot drifters from a summer simulation
dcalcsS = np.load('calcs/coastconn/likelihood/hist-2008-07.npz')
dtracksS = netCDF.Dataset('tracks/2008-07-01T00gc.nc')
xg = dtracksS.variables['xg'][:]; yg = dtracksS.variables['yg'][:]
# xg = xg[::dd,:]; yg = yg[::dd,:]
ind = xg==-1

xpS, ypS, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xpS[ind] = np.nan; ypS[ind] = np.nan

# drifter ids
ids = dcalcsS['ids'].item()
idsuse = []
[idsuse.extend(ids[0,5,i]) for i in xrange(342)]
idsuseS = set(idsuse) # unique only
# set up a set of indices with respect to all drifters with 
# coastal ones as True
indsS = np.zeros(xpS.shape[0]).astype(bool)
indsS[list(idsuseS)] = True

# distance along the coast boxes
dist = np.zeros(len(pathsxy))
verts0 = pathsxy[0].vertices
for i, path in enumerate(pathsxy):
    verts1 = path.vertices
    dist[i:] += np.sqrt((verts1[0,0]-verts0[0,0])**2+(verts1[0,1]-verts0[0,1])**2)
    verts0 = verts1.copy()
dist /= 1000 # convert to km

fig = plt.figure(figsize=(13, 6.6125))
fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
# winter
ax = fig.add_subplot(1,2,1)
ax.set_title('Winter')
# plot background
tracpy.plotting.background(grid, ax=ax, mers=np.arange(-100, -80, 2))
# drifters
ax.plot(xpW[::dd,:].T, ypW[::dd,:].T, '0.4', lw=1, alpha=0.6)
ax.plot(xpW[::dd,0].T, ypW[::dd,0].T, 'o', color='0.4', ms=5)
# drifters that reach the coastline
ax.plot(xpW[indsW,:][::dd,:].T, ypW[indsW,:][::dd,:].T, 'r', lw=2, alpha=0.6)
ax.plot(xpW[indsW,:][::dd,0].T, ypW[indsW,:][::dd,0].T, 'o', color='r', ms=8)
# Plot outline of coastal region
plt.plot(verts_outerxy.vertices[:,0], verts_outerxy.vertices[:,1], 'r-', lw=5, alpha=0.5, zorder=15)
# Plot boxes of coastal region and label
for i, path in enumerate(pathsxy):
    if (i==0):
        ax.text(path.vertices[:,0].min()-30000, path.vertices[:,1].min(), '%.0f' % dist[i], fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif  abs(dist[i] - 250)<2:
        ax.text(path.vertices[:,0].min()-50000, path.vertices[:,1].min()+50000, '%.0f' % 250, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif  abs(dist[i] - 500)<2:
        ax.text(path.vertices[:,0].min()-50000, path.vertices[:,1].min()+50000, '%.0f' % 500, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif abs(dist[i] - 750)<2:
        ax.text(path.vertices[:,0].min()-10000, path.vertices[:,1].min()+50000, '%.0f' % 750, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif abs(dist[i] - 1000)<4:
        ax.text(path.vertices[:,0].min(), path.vertices[:,1].min()+25000, '%.0f' % 1000, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif abs(dist[i] - 1250)<1:
        ax.text(path.vertices[:,0].min(), path.vertices[:,1].min()+25000, '%.0f' % 1250, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))
    elif abs(dist[i] - 1500)<1:
        ax.text(path.vertices[:,0].min()-50000, path.vertices[:,1].min()+50000, '%.0f' % 1500, fontsize=12, 
                color='0.2', bbox=dict(facecolor='w', edgecolor='none', boxstyle='round'))

# summer
ax = fig.add_subplot(1,2,2)
ax.set_title('Summer')
# plot background
tracpy.plotting.background(grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
# drifters
ax.plot(xpS[::dd,:].T, ypS[::dd,:].T, '0.4', lw=1, alpha=0.6)
ax.plot(xpS[::dd,0].T, ypS[::dd,0].T, 'o', color='0.4', ms=5)
# drifters that reach the coastline
ax.plot(xpS[indsS,:][::dd,:].T, ypS[indsS,:][::dd,:].T, 'r', lw=2, alpha=0.6)
ax.plot(xpS[indsS,:][::dd,0].T, ypS[indsS,:][::dd,0].T, 'o', color='r', ms=8)
# Plot outline of coastal region
plt.plot(verts_outerxy.vertices[:,0], verts_outerxy.vertices[:,1], 'r-', lw=5, alpha=0.5, zorder=15)

fig.savefig('figures/coastconn/likelihood/sampledrifters.png', bbox_inches='tight')
