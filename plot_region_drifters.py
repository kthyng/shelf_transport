'''
Plot drifters from winter and summer transport regions, interannual. For paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
from glob import glob
import netCDF4 as netCDF
import op
import pdb
import matplotlib as mpl


shelf_depth = 20
# region = 3
whichseason = 'winter'
cross = False  # to plot the ones that do or do not cross shelf

if whichseason == 'winter':
    region = 3
elif whichseason == 'summer':
    region = 4

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

if cross:
    color = 'r'
else:
    color = 'darkcyan'

grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)

d = np.load('calcs/xyp0.npz')
xp0 = d['xp0']; yp0 = d['yp0']
lonp0, latp0 = grid['basemap'](xp0, yp0, inverse=True)
# X, Y = np.meshgrid(op.resize(d['xe'],0), op.resize(d['ye'],0))
fh = grid['trir'].nn_interpolator(grid['h'].flatten())
depths = fh(xp0, yp0)
# flonr = grid['trir'].nn_interpolator(grid['lonr'].flatten())
# lons = flonr(X,Y)
# flatr = grid['trir'].nn_interpolator(grid['latr'].flatten())
# lats = flatr(X,Y)
d.close()

ishallow = depths < shelf_depth
ideep = depths > shelf_depth

# spacing for drifters
if shelf_depth == 100:
    ds = 500
    ishelf_depth = 2
elif shelf_depth == 50:
    ds = 250
    ishelf_depth = 1
elif shelf_depth == 20:
    ds = 100
    ishelf_depth = 0

years = np.arange(2004, 2015)

# Calculate the transport for each year for the winter wind transport area.
# Want places where the depth is less than 100 meters and west of 90 deg.
if region == 1: # shelf bend
    iwind = ideep * (lonp0<-95) * (lonp0>-97) * (latp0>26) * (latp0<28)
elif region == 2: # Mississippi region
    iwind = ishallow * (lonp0<-88) * (lonp0>-89.3) * (latp0>28.7) * (latp0<30)
elif region == 3: # winter wind transport region
    iwind = ishallow * (lonp0<-90) * (latp0>27.5)
elif region == 4: # summer river transport region
    iwind = ishallow * (lonp0<-90) * (latp0>27.5)
# import pdb; pdb.set_trace()
# iwind = iwind[::ds]

fig, axarr = plt.subplots(4,3)
fig.set_size_inches(8.7, 11.5)
fig.subplots_adjust(left=0.008, bottom=0.1, right=1.0, top=0.98, wspace=0.005, hspace=0.1)

# fig, axarr = plt.subplots(2,4)
# fig.set_size_inches(13.4, 6.6125)
# fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

for i, ax in enumerate(axarr.flatten()):

    if i != 11:
        if whichseason == 'winter':
            Files = glob('tracks/' + str(years[i]) + '-0[1-2]-??T0?gc.nc')
        elif whichseason == 'summer':
            Files = glob('tracks/' + str(years[i]) + '-0[7-8]-??T0?gc.nc')

    # Titles for subplots
    if i==10:#4:
        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
            pars=np.arange(20, 36, 2), outline=True, parslabels=[0, 1, 0, 0])
        # ax.set_title(str(2004+i))
    elif i==11:#7:
        # ax.set_frame_on(False)
        ax.set_axis_off()
    else:
        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
            pars=np.arange(20, 36, 2), outline=True,
            merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
        # ax.set_title(str(2004+i))

    if i!=11:
        ax.text(0.07, 0.88, str(2004+i), transform=ax.transAxes)

    ax.set_frame_on(False)

    if i == 11:
        continue

    # loop through tracks and add trajectories from region to plot
    for File in Files[::6]:
        print(File)
        d = netCDF.Dataset(File)
        xg = d.variables['xg'][:]; yg = d.variables['yg'][:]
        # import pdb; pdb.set_trace()
        nanind = (xg[iwind,:][::ds,:]==-1)
        xp, yp, _ = tracpy.tools.interpolate2d(xg[iwind,:][::ds,:], yg[iwind,:][::ds,:], grid, 'm_ij2xy')
        xp[nanind] = np.nan; yp[nanind] = np.nan

        # indices of those that do or do not cross shelf, as desired
        shelffile = 'calcs/shelfconn/' + File.split('/')[-1][:-3] + '.npz'  # file of cross-shelf calculations
        # indices of drifters that cross 100 m isobath
        if cross:
            icross = ~np.isnan(np.load(shelffile)['cross'][ishelf_depth, iwind][::ds])
        else:
            icross = np.isnan(np.load(shelffile)['cross'][ishelf_depth, iwind][::ds])

        ax.plot(xp[icross,:].T, yp[icross,:].T, color=color, alpha=0.2, lw=0.1)

    # pdb.set_trace()

if cross:
    plt.savefig('figures/plot_region_' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-cross_lowres.png', bbox_inches='tight', dpi=100)
    plt.savefig('figures/plot_region_' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-cross.png', bbox_inches='tight', dpi=300)
else:
    plt.savefig('figures/plot_region_' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-notcross_lowres.png', bbox_inches='tight', dpi=100)
    plt.savefig('figures/plot_region_' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-notcross.png', bbox_inches='tight', dpi=300)
