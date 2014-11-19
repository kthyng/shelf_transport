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
import op
from matplotlib import ticker
from matplotlib.mlab import find


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

whichtime = 'interannual' # 'seasonal' or 'interannual'
whicharea = 'winter' # 'winter' or 'summer'

# which drifters to plot for summer and winter, respectively
inds = np.load('calcs/summer-drifter-indices.npz')['inds']
indw = np.load('calcs/winter-drifter-indices.npz')['inds']

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# decimate the drifter indices
ddi = 1000

if plottracks:

    if whichtime == 'seasonal':

        if whicharea == 'winter':

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

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area.png', bbox_inches='tight')


        elif whicharea == 'summer':

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

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area.png', bbox_inches='tight')


    elif whichtime == 'interannual':

        if whicharea == 'winter':

            # WINTER TRANSPORT AREA, showing winter interannual

            fig, axarr = plt.subplots(2,4)
            fig.set_size_inches(13.4, 6.6125)
            fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

            for i, ax in enumerate(axarr.flatten()):

                yr = 2004+i
                Files = glob.glob('tracks/' + str(yr) + '-0[1,2]-*gc.nc')

                # Titles for subplots
                if i==4:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                elif i==7:
                    ax.set_frame_on(False)
                    ax.set_axis_off()
                    continue
                else:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                        merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])

                ax.set_title(str(yr))
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

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area.png', bbox_inches='tight')
        
        elif whicharea == 'summer':

            # SUMMER TRANSPORT AREA, showing summer interannual

            fig, axarr = plt.subplots(2,4)
            fig.set_size_inches(13.4, 6.6125)
            fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

            for i, ax in enumerate(axarr.flatten()):

                yr = 2004+i
                
                # Titles for subplots
                if i==4:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                elif i==7:
                    ax.set_frame_on(False)
                    ax.set_axis_off()
                    continue
                else:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                        merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])

                ax.set_title(str(yr))
                Files = glob.glob('tracks/' + str(yr) + '-0[7,8]-*gc.nc')
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

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area.png', bbox_inches='tight')
        

elif plothist:

    if whichtime == 'seasonal':

        cmap = 'YlGn'

        shelf_depth = 100

        # Calculate xrange and yrange for histograms
        XPrange = [grid['xpsi'].min(), grid['xpsi'].max()]
        YPrange = [grid['ypsi'].min(), grid['ypsi'].max()]
        # # Convert projected range into grid space to save conversion
        # XGrange, YGrange, _ = tracpy.tools.interpolate2d(XPrange, YPrange, grid, 'd_xy2ij')

        # Number of bins to use in histogram
        bins = (100,100) #(30,30)

        H = np.zeros(bins) # initialize

        if whicharea == 'winter':

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

                for i,File in enumerate(Files):
                    # print File
                    d = netCDF.Dataset(File)
                    # pdb.set_trace()
                    # reading in a subset of indices from the beginning is prohibitively slow
                    xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
                    yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
                    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
                    nind = xg==-1
                    del(xg,yg)
                    # xg[nind] = np.nan; yg[nind] = np.nan
                    xp[nind] = np.nan; yp[nind] = np.nan
                    d.close()

                    # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                    # Htemp, xe, ye = np.histogram2d(xg.flatten(), yg.flatten(), bins=bins, range=[[XGrange[0], XGrange[1]], [YGrange[0], YGrange[1]]])
                    Htemp, xe, ye = np.histogram2d(xp.flatten(), yp.flatten(), bins=bins, range=[[XPrange[0], XPrange[1]], [YPrange[0], YPrange[1]]])
                    H = np.nansum( np.vstack((H[np.newaxis,:,:], Htemp[np.newaxis,:,:])), axis=0)

                XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))
                hist, bin_edges = np.histogram(H.flat, bins=100) # find # of occurrences of histogram bin values
                n = np.cumsum(hist)
                Hmax = bin_edges[find(n<(n.max()-n.min())*.8+n.min())[-1]] # take the 80% of histogram occurrences as the max instead of actual max since too high
                locator = ticker.MaxNLocator(11)
                locator.create_dummy_axis()
                locator.set_bounds(0, 1) 
                levels = locator()
                extend = 'max'
                H = H/Hmax
                mappable = ax.contourf(XE, YE, H.T, cmap=cmap, levels=levels, extend=extend)
                ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

                # outline the area where drifters started
                d = np.load('calcs/winter-contour-pts.npz')
                ax.plot(d['x'], d['y'], 'k', lw=3)
                d.close()

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area-hist.png', bbox_inches='tight')

        elif whicharea == 'summer':

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

                for i,File in enumerate(Files):
                    # print File
                    d = netCDF.Dataset(File)
                    # pdb.set_trace()
                    # reading in a subset of indices from the beginning is prohibitively slow
                    xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
                    yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
                    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
                    nind = xg==-1
                    del(xg,yg)
                    # xg[nind] = np.nan; yg[nind] = np.nan
                    xp[nind] = np.nan; yp[nind] = np.nan
                    d.close()

                    # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                    # Htemp, xe, ye = np.histogram2d(xg.flatten(), yg.flatten(), bins=bins, range=[[XGrange[0], XGrange[1]], [YGrange[0], YGrange[1]]])
                    Htemp, xe, ye = np.histogram2d(xp.flatten(), yp.flatten(), bins=bins, range=[[XPrange[0], XPrange[1]], [YPrange[0], YPrange[1]]])
                    H = np.nansum( np.vstack((H[np.newaxis,:,:], Htemp[np.newaxis,:,:])), axis=0)

                XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))
                hist, bin_edges = np.histogram(H.flat, bins=100) # find # of occurrences of histogram bin values
                n = np.cumsum(hist)
                Hmax = bin_edges[find(n<(n.max()-n.min())*.8+n.min())[-1]] # take the 80% of histogram occurrences as the max instead of actual max since too high
                locator = ticker.MaxNLocator(11)
                locator.create_dummy_axis()
                locator.set_bounds(0, 1) 
                levels = locator()
                extend = 'max'
                H = H/Hmax
                mappable = ax.contourf(XE, YE, H.T, cmap=cmap, levels=levels, extend=extend)
                ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

                # outline the area where drifters started
                d = np.load('calcs/summer-contour-pts.npz')
                ax.plot(d['x'], d['y'], 'k', lw=3)
                d.close()

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area-hist.png', bbox_inches='tight')

    elif whichtime == 'interannual':

        cmap = 'YlGn'

        shelf_depth = 100

        # Calculate xrange and yrange for histograms
        XPrange = [grid['xpsi'].min(), grid['xpsi'].max()]
        YPrange = [grid['ypsi'].min(), grid['ypsi'].max()]
        # # Convert projected range into grid space to save conversion
        # XGrange, YGrange, _ = tracpy.tools.interpolate2d(XPrange, YPrange, grid, 'd_xy2ij')

        # Number of bins to use in histogram
        bins = (100,100) #(30,30)

        H = np.zeros(bins) # initialize

        if whicharea == 'winter':

            # WINTER TRANSPORT AREA, showing both winter and summer seasons

            fig, axarr = plt.subplots(2,4)
            fig.set_size_inches(13.4, 6.6125)
            fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

            for i, ax in enumerate(axarr.flatten()):

                # Titles for subplots
                if i==4:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                    ax.set_title(str(2004+i))
                    Files = glob.glob('tracks/' + str(yr) + '-0[1,2]-*gc.nc')
                elif i==7:
                    ax.set_frame_on(False)
                    ax.set_axis_off()
                    continue
                else:
                    yr = 2004+i
                    Files = glob.glob('tracks/' + str(yr) + '-0[1,2]-*gc.nc')
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                        merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    ax.set_title(str(yr))

                ind = indw

                for i,File in enumerate(Files):
                    # print File
                    d = netCDF.Dataset(File)
                    # pdb.set_trace()
                    # reading in a subset of indices from the beginning is prohibitively slow
                    xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
                    yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
                    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
                    nind = xg==-1
                    del(xg,yg)
                    # xg[nind] = np.nan; yg[nind] = np.nan
                    xp[nind] = np.nan; yp[nind] = np.nan
                    d.close()

                    # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                    # Htemp, xe, ye = np.histogram2d(xg.flatten(), yg.flatten(), bins=bins, range=[[XGrange[0], XGrange[1]], [YGrange[0], YGrange[1]]])
                    Htemp, xe, ye = np.histogram2d(xp.flatten(), yp.flatten(), bins=bins, range=[[XPrange[0], XPrange[1]], [YPrange[0], YPrange[1]]])
                    H = np.nansum( np.vstack((H[np.newaxis,:,:], Htemp[np.newaxis,:,:])), axis=0)

                XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))
                hist, bin_edges = np.histogram(H.flat, bins=100) # find # of occurrences of histogram bin values
                n = np.cumsum(hist)
                Hmax = bin_edges[find(n<(n.max()-n.min())*.8+n.min())[-1]] # take the 80% of histogram occurrences as the max instead of actual max since too high
                locator = ticker.MaxNLocator(11)
                locator.create_dummy_axis()
                locator.set_bounds(0, 1) 
                levels = locator()
                extend = 'max'
                H = H/Hmax
                mappable = ax.contourf(XE, YE, H.T, cmap=cmap, levels=levels, extend=extend)
                ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

                # outline the area where drifters started
                d = np.load('calcs/winter-contour-pts.npz')
                ax.plot(d['x'], d['y'], 'k', lw=3)
                d.close()

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area-hist.png', bbox_inches='tight')

        elif whicharea == 'summer':

            # SUMMER TRANSPORT AREA, showing both winter and summer seasons

            fig, axarr = plt.subplots(2,4)
            fig.set_size_inches(13.4, 6.6125)
            fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

            for i, ax in enumerate(axarr.flatten()):

                # Titles for subplots
                if i==4:
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                    ax.set_title(str(2004+i))
                    Files = glob.glob('tracks/' + str(yr) + '-0[7,8]-*gc.nc')
                elif i==7:
                    ax.set_frame_on(False)
                    ax.set_axis_off()
                    continue
                else:
                    yr = 2004+i
                    Files = glob.glob('tracks/' + str(yr) + '-0[7,8]-*gc.nc')
                    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                        merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    ax.set_title(str(yr))

                ind = inds

                for i,File in enumerate(Files):
                    # print File
                    d = netCDF.Dataset(File)
                    # pdb.set_trace()
                    # reading in a subset of indices from the beginning is prohibitively slow
                    xg = d.variables['xg'][:]; xg = xg[ind[::ddi],:]
                    yg = d.variables['yg'][:]; yg = yg[ind[::ddi],:]
                    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
                    nind = xg==-1
                    del(xg,yg)
                    # xg[nind] = np.nan; yg[nind] = np.nan
                    xp[nind] = np.nan; yp[nind] = np.nan
                    d.close()

                    # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                    # Htemp, xe, ye = np.histogram2d(xg.flatten(), yg.flatten(), bins=bins, range=[[XGrange[0], XGrange[1]], [YGrange[0], YGrange[1]]])
                    Htemp, xe, ye = np.histogram2d(xp.flatten(), yp.flatten(), bins=bins, range=[[XPrange[0], XPrange[1]], [YPrange[0], YPrange[1]]])
                    H = np.nansum( np.vstack((H[np.newaxis,:,:], Htemp[np.newaxis,:,:])), axis=0)

                XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))
                hist, bin_edges = np.histogram(H.flat, bins=100) # find # of occurrences of histogram bin values
                n = np.cumsum(hist)
                Hmax = bin_edges[find(n<(n.max()-n.min())*.8+n.min())[-1]] # take the 80% of histogram occurrences as the max instead of actual max since too high
                locator = ticker.MaxNLocator(11)
                locator.create_dummy_axis()
                locator.set_bounds(0, 1) 
                levels = locator()
                extend = 'max'
                H = H/Hmax
                mappable = ax.contourf(XE, YE, H.T, cmap=cmap, levels=levels, extend=extend)
                ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

                # outline the area where drifters started
                d = np.load('calcs/summer-contour-pts.npz')
                ax.plot(d['x'], d['y'], 'k', lw=3)
                d.close()

            fig.savefig('figures/drifters/transport-areas/' + whichtime + '-' + whicharea + '-area-hist.png', bbox_inches='tight')
