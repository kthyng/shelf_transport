'''
Plot a set of drifters in map view, in time.
'''

import matplotlib.pyplot as plt
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from matplotlib.mlab import find
import pdb
import numpy as np
import matplotlib as mpl
import os

mpl.rcParams.update({'font.size': 20})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

# read in grid
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# loc = '/home/kthyng/shelf/grid.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# whether to do tails on drifters or not (don't with low decimation)
dotails = False # True or False

# Read in drifter tracks
dd = 1 # 500 # drifter decimation
startdate = '2004-07-15T00'
d = netCDF.Dataset('tracks/' + startdate + 'gc.nc')
xg = d.variables['xg'][::dd,:]
yg = d.variables['yg'][::dd,:]
ind = (xg == -1)
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xp[ind] = np.nan; yp[ind] = np.nan
tp = d.variables['tp'][:]
# d.close()

# Plot drifters, starting 5 days into simulation
# 2 days for one tail, 3 days for other tail
# t = tp-tp[0]
days = (tp-tp[0])/(3600.*24)
dates = netCDF.num2date(tp, d.variables['tp'].units)
# Find indices relative to present time
i5daysago = 0 # keeps track of index 5 days ago
i2daysago = find(days>=3)[0] # index for 2 days ago, which starts as 3 days in
if dotails:
    i5days = find(days>=5)[0] # index for 5 days in
else:
    i5days = 0 # start at the beginning
nt = tp.size # total number of time indices
# for i in np.arange(0,nt+1,5):

dirname = 'figures/drifters/dd' + str(dd) + '/' + startdate
if not os.path.exists(dirname):
    os.makedirs(dirname)

for i in np.arange(i5days,nt+1,5):

    fname = dirname + '/' + dates[i].isoformat()[:-6] + '.png'

    if os.path.exists(fname):
        # Update indices
        i5daysago += 5
        i2daysago += 5
        continue

    # Plot background
    fig = plt.figure(figsize=(18,14))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax)

    if dotails:
        # Plot 5 days ago to 2 days ago
        ax.plot(xp[:,i5daysago:i2daysago].T, yp[:,i5daysago:i2daysago].T, color='0.6', lw=2)

        # Plot 0-2 day tail
        ax.plot(xp[:,i2daysago:i].T, yp[:,i2daysago:i].T, color='0.3', lw=3)

        # Plot drifter locations
        ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='r', ms=10)

    else:

        # Plot drifter locations
        ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='g', ms=2, mec='None')

    # Time
    ax.text(0.075, 0.95, dates[i].isoformat()[:-6], transform=ax.transAxes, fontsize=20)

    # Drifter legend
    if dotails:
        ax.plot(0.0895, 0.9, 'or', ms=10, transform=ax.transAxes) # drifter head
        ax.plot([0.075, 0.1], [0.875, 0.875], '0.3', lw=3, transform=ax.transAxes) # drifter tail #1
        ax.plot([0.075, 0.1], [0.85, 0.85], '0.5', lw=2, transform=ax.transAxes) # drifter tail #2
        ax.text(0.125, 0.89, 'Drifter location', color='r', transform=ax.transAxes, fontsize=16)
        ax.text(0.125, 0.866, '2 days prior', color='0.3', transform=ax.transAxes, fontsize=16)
        ax.text(0.125, 0.842, '5 days prior', color='0.5', transform=ax.transAxes, fontsize=16)
    else:
        ax.plot(0.0895, 0.9, 'o', color='g', ms=10, transform=ax.transAxes) # drifter head
        ax.text(0.125, 0.89, 'Drifter location', color='g', transform=ax.transAxes, fontsize=16)

    # Update indices
    i5daysago += 5
    i2daysago += 5

    # pdb.set_trace()

    fig.savefig(fname, bbox_inches='tight')

    plt.close()

# plt.show()
