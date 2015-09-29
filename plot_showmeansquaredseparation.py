'''
Demonstrate mean squared separation distance with plot.
'''

import numpy as np
import tracpy
import tracpy.plotting
import tracpy.calcs
import netCDF4 as netCDF
import matplotlib.pyplot as plt
import matplotlib as mpl


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


loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid_orig = tracpy.inout.readgrid(loc, usebasemap=True)
grid_small = tracpy.inout.readgrid(loc, usebasemap=True, urcrnrlon=-95, urcrnrlat=27.5,
                                   llcrnrlon=-98, llcrnrlat=24.5)

d = netCDF.Dataset('tracks/2004-01-01T00gc.nc')
xg = d.variables['xg'][100000:100050, :]
yg = d.variables['yg'][100000:100050, :]
# convert to lat/lon
lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid_orig, 'm_ij2ll')
# convert to project coords in smaller grid for plotting
xp, yp = grid_small['basemap'](lonp, latp)
t = d.variables['tp'][:]
days = (t - t[0])/(3600.*24)

# Calculate relative dispersion between pairs
D2p1, _, _ = tracpy.calcs.rel_dispersion(xp[3:5,:], yp[3:5,:], r=[0,2], squared=True, spherical=False)
D2p2, _, _ = tracpy.calcs.rel_dispersion(xp[4:6,:], yp[4:6,:], r=[0,2], squared=True, spherical=False)
D2both, _, _ = tracpy.calcs.rel_dispersion(xp[3:6,:], yp[3:6,:], r=[0,2], squared=True, spherical=False)

# Make a figure with the drifter pairs and then relative dispersion above
fig = plt.figure(figsize=(12,9))
ax1 = plt.subplot2grid((4, 2), (0, 0), colspan=2)  # relative dispersion
ax2 = plt.subplot2grid((4, 2), (1, 0), rowspan=3)  # drifter pair #1
ax3 = plt.subplot2grid((4, 2), (1, 1), rowspan=3)  # drifter pair #2
fig.subplots_adjust(left=0.04, bottom=0.06, right=0.98, top=0.94, wspace=0.05, hspace=0.29)

# plot mean squared separation distance
for i in np.arange(0, days.size, 10):
    ax1.set_xlim(0, 30)
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax1.plot(days[i], D2p1[i], 'bo', ms=10)
    ax1.plot(days[i], D2p2[i], 'go', ms=10)
    ax1.plot(days[i], D2both[i], 'o', color='0.3', ms=10)
    ax1.plot(days[:i+1], D2p1[:i+1], 'b', days[:i+1], D2p2[:i+1], 'g', days[:i+1], D2both[:i+1], '0.3', lw=5)
    if i < 385:  # use smaller y axis
        ax1.set_ylim(0, 25000)
    else:  # use full y axis
        ax1.set_ylim(0, 90000)
    ax1.set_title('Mean squared separation distance [km$^2\!$]')
    ax1.set_xlabel('Time [days]')

    # plot drifters that stay close together
    tracpy.plotting.background(grid_small, fig=fig, ax=ax2, outline=[0, 0, 0, 0])
    ax2.plot(xp[3, 0].T, yp[3, 0].T, 'bo', ms=15, alpha=0.4)
    ax2.plot(xp[3:5, i].T, yp[3:5, i].T, 'bo', ms=15, alpha=1.0)
    ax2.plot(xp[3:5, :i].T, yp[3:5, :i].T, 'b-', lw=4)

    # plot drifters that separate
    tracpy.plotting.background(grid_small, fig=fig, ax=ax3, outline=[0, 0, 0, 0], parslabels=[0, 0, 0, 0],
                               mers=np.arange(-97, -94))
    ax3.plot(xp[4, 0].T, yp[4, 0].T, 'go', ms=15, alpha=0.4)
    ax3.plot(xp[4:6, i].T, yp[4:6, i].T, 'go', ms=15, alpha=1.0)
    ax3.plot(xp[4:6, :i].T, yp[4:6, :i].T, 'g-', lw=4)

    fig.savefig('figures/showmeansquaredseparation/' + str(i).zfill(3) + '.png', bbox_inches='tight', dpi=80)

    ax1.clear()
    ax2.clear()
    ax3.clear()
