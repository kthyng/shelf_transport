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
import tracpy.calcs

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
startdate = '2010-01-15T00'
d = netCDF.Dataset('tracks/' + startdate + 'gc.nc')
xg = d.variables['xg'][::dd,:]
yg = d.variables['yg'][::dd,:]
ind = (xg == -1)
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xp[ind] = np.nan; yp[ind] = np.nan
tp = d.variables['tp'][:]
# d.close()

# txla output
nc = netCDF.Dataset(loc)
datestxla = netCDF.num2date(nc.variables['ocean_time'][:], nc.variables['ocean_time'].units)

if not dotails:
    
    # Find indices of drifters by starting depth
    depthp = tracpy.calcs.Var(xg[:,0], yg[:,0], tp, 'h', nc) # starting depths of drifters
    # near-shore
    # ind10 = depthp<=10
    ind20 = depthp<=20
    ind50 = (depthp>20)*(depthp<=50)
    ind100 = (depthp>50)*(depthp<=100)
    # offshore
    ind500 = (depthp>100)*(depthp<=500)
    ind3500 = depthp>500

    # colors for drifters
    rgb = plt.cm.get_cmap('winter_r')(np.linspace(0,1,6))[:-1,:3] # skip last entry where it levels off in lightness

    # to plot colorbar
    gradient = np.linspace(0, 1, 6)[:-1]
    gradient = np.vstack((gradient, gradient))

    ms = 1.5 # markersize


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

    if not dotails:
        itxla = np.where(datestxla==dates[i])[0][0] # find time index to use for model output
        salt = nc.variables['salt'][itxla,-1,:,:] # surface salinity

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
        # ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='g', ms=2, mec='None')
        ax.plot(xp[ind20,i].T, yp[ind20,i].T, 'o', color=rgb[0,:], ms=ms, mec='None')
        ax.plot(xp[ind50,i].T, yp[ind50,i].T, 'o', color=rgb[1,:], ms=ms, mec='None')
        ax.plot(xp[ind100,i].T, yp[ind100,i].T, 'o', color=rgb[2,:], ms=ms, mec='None')
        ax.plot(xp[ind500,i].T, yp[ind500,i].T, 'o', color=rgb[3,:], ms=ms, mec='None')
        ax.plot(xp[ind3500,i].T, yp[ind3500,i].T, 'o', color=rgb[4,:], ms=ms, mec='None')

        # Overlay surface salinity
        ax.contour(grid['xr'].T, grid['yr'].T, salt, [33], colors='0.1', zorder=12, linewidths=2)

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
        cax = fig.add_axes([0.2, 0.8, 0.15, 0.02])
        cax.imshow(gradient, aspect='auto', interpolation='none', cmap=plt.get_cmap('winter_r'))
        cax.tick_params(axis='y', labelleft=False, left=False, right=False)
        cax.tick_params(axis='x', top=False, bottom=False, labelsize=15)
        cax.set_xticks(np.arange(-0.5, 5, 1.0))
        cax.set_xticklabels(('0', '20', '50', '100', '500', '3500'))
        cax.set_title('Initial drifter depth [m]', fontsize=16)
        # legend for contour
        ax.plot([0.075, 0.1], [0.81, 0.81], '0.1', lw=2, transform=ax.transAxes)
        ax.text(0.123, 0.802, '33 salinity contour', color='0.1', transform=ax.transAxes, fontsize=16)


    # Update indices
    i5daysago += 5
    i2daysago += 5

    fig.savefig(fname, bbox_inches='tight')

    plt.close()
