'''
Plot some drifter tracks and location histograms, as determined by input indices.
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
import matplotlib as mpl
import pdb


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


plottracks = True # plot drifter trajectories
plothist = False # plot drifter locations at a time as a histogram

# which drifters to plot for summer and winter, respectively
inds = np.load('calcs/summer-drifter-indices.npz')['inds']
indw = np.load('calcs/winter-drifter-indices.npz')['inds']

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# decimate the drifter indices
ddi = 1000

if plottracks:

    # WINTER TRANSPORT AREA, showing both winter and summer seasons

    fig, axarr = plt.subplots(1,2)
    fig.set_size_inches(13.675, 6.6125)
    fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)

    for i, ax in enumerate(axarr):

        # Titles for subplots
        if i==0:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title('Winter')
            Files = glob.glob('tracks/20??-0[1,2]-*gc.nc')
            ind = indw 
        elif i==1:
            tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            Files = glob.glob('tracks/20??-0[7,8]-*gc.nc')
            ind = indw


        xpf = []; ypf = []
        for i,File in enumerate(Files):
            # print File
            d = netCDF.Dataset(File)
            # pdb.set_trace()
            # reading in a subset of indices from the beginning is prohibitively slow
            xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
            yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
            nind = xg==-1
            xp[nind] = np.nan; yp[nind] = np.nan
            d.close()

            ax.plot(xp.T, yp.T, '0.3', lw=0.5, alpha=0.5) # tracks
            xpft, ypft = tracpy.tools.find_final(xp, yp)
            xpf.append(xpft); ypf.append(ypft)

        # plot the ending locations at the end so we can see them
        ax.plot(xpf, ypf, 'ro', ms=3)
        # outline the area where drifters started
        d = np.load('calcs/winter-contour-pts.npz')
        ax.plot(d['x'], d['y'], 'k', lw=3)
        d.close()

    fig.savefig('figures/drifters/transport-areas/winter-area.png', bbox_inches='tight')



    # SUMMER TRANSPORT AREA, showing both winter and summer seasons

    fig, axarr = plt.subplots(1,2)
    fig.set_size_inches(13.675, 6.6125)
    fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)

    for i, ax in enumerate(axarr):

        # Titles for subplots
        if i==0:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title('Winter')
            Files = glob.glob('tracks/20??-0[1,2]-*gc.nc')
            ind = inds 
        elif i==1:
            tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            Files = glob.glob('tracks/20??-0[7,8]-*gc.nc')
            ind = inds


        xpf = []; ypf = []
        for i,File in enumerate(Files):
            # print File
            d = netCDF.Dataset(File)
            # pdb.set_trace()
            # reading in a subset of indices from the beginning is prohibitively slow
            xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
            yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
            nind = xg==-1
            xp[nind] = np.nan; yp[nind] = np.nan
            d.close()

            ax.plot(xp.T, yp.T, '0.3', lw=0.5, alpha=0.5) # tracks
            xpft, ypft = tracpy.tools.find_final(xp, yp)
            xpf.append(xpft); ypf.append(ypft)

        # plot the ending locations at the end so we can see them
        ax.plot(xpf, ypf, 'ro', ms=3)
        # outline the area where drifters started
        d = np.load('calcs/summer-contour-pts.npz')
        ax.plot(d['x'], d['y'], 'k', lw=3)
        d.close()

    fig.savefig('figures/drifters/transport-areas/summer-area.png', bbox_inches='tight')

        # plot_colorbar(fig, mappable, 'diff', ticks=ticks)
