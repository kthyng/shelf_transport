'''
Present the likelihood that oil spilled in an area will reach the coastline.
Use analysis already calculated in find_coastal_path_connectivity.
Also look at the vulnerability of the coastline for the boxes.
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
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar
import matplotlib.patches as Patches
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from collections import defaultdict
import cmocean.cm as cmo

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


# grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
# grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)
# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
proj = tracpy.tools.make_proj('nwgom')
grid = tracpy.inout.readgrid('../grid.nc', proj)

# load in initial drifter starting locations in grid space
d = np.load('calcs/xyp0.npz')
xp0 = d['xp0']; yp0 = d['yp0']
d.close()

bins = (100,100)
# Calculate xrange and yrange for histograms
Xrange = [grid.x_psi.min(), grid.x_psi.max()]
Yrange = [grid.y_psi.min(), grid.y_psi.max()]

# Save a histogram of the number of drifters that started in each bin
Hstartfile = 'calcs/coastconn/likelihood/Hstart.npz'
if not os.path.exists(Hstartfile):
    Hstart, xe, ye = np.histogram2d(xp0, yp0, bins=bins,
                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
    np.savez(Hstartfile, Hstart=Hstart, xe=xe, ye=ye)
else:
    d = np.load(Hstartfile)
    Hstart = d['Hstart']
    xe = d['xe']; ye = d['ye']
    d.close()


# # Find indices of all drifters that start in the coastal boxes
# # start indices of drifters that start in each box path
pts = np.load('calcs/alongcoastconn/inds-in-coast-paths.npz')['pts']
# pt = [] # aggregated indices of drifters that start in coast boxes
# [pt.extend(pts[j]) for j in xrange(len(pts))]
# xp0coast = xp0[pt]; yp0coast = yp0[pt]

# # Get and read in shipping lanes data
# # get data
# # info: https://www.data.boem.gov/homepg/pubinfo/repcat/arcinfo/zipped/gomr_fairways.htm
# url = 'http://www.data.boem.gov/homepg/pubinfo/repcat/arcinfo/zipped/fairway.zip'
# dirname = 'fairway'
# if not os.path.exists(dirname):
#     os.mkdir('fairway')
# if not os.path.exists('fairway.zip'):
#     os.system('wget ' + url)
# os.system('unzip -d fairway fairway.zip')

# read in shipping lanes data
proj.readshapefile('fairway/fairway', 'fairway')


def plot_interannual():
    '''
    Plot interannual comparison of likelihood, either overall or just certain parts.
    '''

    cmap = 'YlGn'
    log = False
    season = 'summer' # 'winter' or 'summer'
    # zoomed = True
    whichboxes = 'both' # 'all', 'porta', 'galveston', 'both'
    # which boxes along coast to use for vulnerability. 0:342 for all
    # Port A:
    if whichboxes=='all':
        boxes = np.arange(0,342)
        whichH = 'Hall' # Hall for histogram for entire coastline at once; H for by coast box
        zoomed = False
    elif whichboxes=='all2':
        boxes = np.arange(0,342)
        whichH = 'H' # use coast box histograms instead of combined
        zoomed = False
    elif whichboxes=='porta':
        boxes = np.arange(103,123)
        whichH = 'H'
        zoomed = True
        # limits for zoomed box
        if season == 'winter':
        # x1, x2, y1, y2 = 86000, 340800, 465000, 715000
            x1, x2, y1, y2 = 86000, 277100, 527500, 715000
        elif season == 'summer':
            x1, x2, y1, y2 = 86000, 277100, 465000, 652500
        zoom = 2.5
        plume = False
    elif whichboxes=='galveston':
        boxes = np.arange(160,180)
        whichH = 'H'
        zoomed = True
        x1, x2, y1, y2 = 277100, 531900, 660000, 810000
        zoom = 2.75
        plume = False
        # x1, x2, y1, y2 = 277100, 531900, 560000, 810000
        # zoom = 2.0
    elif whichboxes=='both':
        boxes = np.arange(103,180)
        whichH = 'H'
        zoomed = True
        x1, x2, y1, y2 = 86000, 531900, 465000, 810000
        zoom = 1.5

    d = np.load('calcs/coastpaths.npz')
    pathsxy = d['pathsxy']
    d.close()

    ## Read in files ##
    filename = 'calcs/coastconn/likelihood/hist-' + season + 'interannual-' + whichboxes + '.npz'
    if not os.path.exists(filename):
        base = 'calcs/coastconn/likelihood/hist-'
        Files = []
        if season == 'winter':
            months = '[1-2]'
        elif season == 'summer':
            months = '[7-8]'
        Files.append(glob.glob(base + '2004-0' + months + '.npz'))
        Files.append(glob.glob(base + '2005-0' + months + '.npz'))
        Files.append(glob.glob(base + '2006-0' + months + '.npz'))
        Files.append(glob.glob(base + '2007-0' + months + '.npz'))
        Files.append(glob.glob(base + '2008-0' + months + '.npz'))
        Files.append(glob.glob(base + '2009-0' + months + '.npz'))
        Files.append(glob.glob(base + '2010-0' + months + '.npz'))
        Files.append(glob.glob(base + '2011-0' + months + '.npz'))
        Files.append(glob.glob(base + '2012-0' + months + '.npz'))
        Files.append(glob.glob(base + '2013-0' + months + '.npz'))
        Files.append(glob.glob(base + '2014-0' + months + '.npz'))
        H = np.zeros((len(Files),6,100,100))
        ndbox = np.zeros((len(Files),6,342))
        for i,files in enumerate(Files): # winter and summer
            numfiles = 0
            for File in files: # months/years within winter or summer
                print File
                d = np.load(File)
                numfiles += d['numfiles']
                idstemp = d['ids'].item()
                for k in xrange(d['numfiles']): # loop over simulations in file which don't want to mix up
                    for j in xrange(6): #advection days
                        # add together the drifters for the desired boxes
                        ids = []
                        [ids.extend(idstemp[k,j,box]) for box in boxes]
                        ids = list(set(ids)) # eliminate nonunique drifters
                        Htemp, _, _ = np.histogram2d(xp0[ids], yp0[ids], bins=bins,
                                        range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                        H[i,j] += Htemp
                days = d['days']
                xe = d['xe']; ye = d['ye']
                ndbox[i,:,:] += d['ndbox']
                d.close()
            # Divide by number of starting drifters
            H[i,:,:,:] /= (numfiles*Hstart)
        np.savez(filename, H=H, days=days, xe=xe, ye=ye, ndbox=ndbox)

    else:
        d = np.load(filename)
        H = d['H']; days = d['days']; xe = d['xe']; ye = d['ye']; ndbox = d['ndbox']
        d.close()
    ####

    # limit boxes
    pathsxy = pathsxy[boxes]
    ndbox = ndbox[:,:,boxes]

    XE, YE = np.meshgrid(tracpy.op.resize(xe, 0), tracpy.op.resize(ye, 0))

    # Loop through advection days
    # for j in xrange(H.shape[1]):
    j = H.shape[1]-1

    ## Plot setup ##
    fig, axarr = plt.subplots(4,3)
    fig.set_size_inches(8.7, 11.5)
    fig.subplots_adjust(left=0.008, bottom=0.1, right=1.0, top=0.98, wspace=0.005, hspace=0.1)

    for i, ax in enumerate(axarr.flatten()):
       # Titles for subplots
        if i==10:#4:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3),
                pars=np.arange(20, 36, 2), outline=False, parslabels=[0, 1, 0, 0])
        elif i==11:#7:
            ax.set_axis_off()
        else:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3),
                pars=np.arange(20, 36, 2), outline=False,
                merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
        if i!=11:
            # if whichH=='Hall':
            levels = np.linspace(0,100,11)
            var = H[i,j].T*100.
            # elif whichH=='H':
            #     levels = np.linspace(0,1,11)
            #     var = H[i,j].T/H[:,j].max()
            ax.text(0.07, 0.88, str(2004+i), transform=ax.transAxes)
            if log:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels, norm=colors.LogNorm())
            else:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels)

                if plume: # add isohalines from surface salinity
                    dsalt = np.load('calcs/coastconn/likelihood/salt' + str(2004+i) + '.npz')
                    salt = dsalt['salt']; dates = dsalt['dates']
                    dsalt.close()

                    # pdb.set_trace()

                    # for k in xrange(salt.shape[0]):
                    #     # pdb.set_trace()
                    #     ax.contour(grid['xr'], grid['yr'], salt[k,:,:].T, [32], colors='b', linewidths=0.05, alpha=0.3)
                    ax.contour(grid.xr, grid.yr, salt.mean(axis=0).T, [32], colors='b', linewidths=1, alpha=0.8)

                # Add on vulnerability of coastline
                # Need to plot the coast boxes as patches and color them according to vulnerability level
                # http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
                # we need to normalize the data to 0..1 for the full
                # range of the colormap
                fracs = ndbox[i,j,:].astype(float)/ndbox[:,j,:].max() # max across years
                norm = colors.Normalize(fracs.min(), fracs.max())

                # Save patches together
                patches = []
                for path in pathsxy:
                    patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5))

                # assign shades of colormap to the patches according to values, and plot
                for thisfrac, thispatch in zip(fracs, patches):
                    color = cmo.matter(norm(thisfrac))
                    thispatch.set_facecolor(color)
                    ax.add_patch(thispatch)

                if lanes:
                    # plot shipping lanes
                    for entry in proj.fairway:
                        entry = np.asarray(entry)
                        ax.plot(entry[:, 0], entry[:, 1], 'r')

                if zoomed:# and j==H.shape[1]-1: # magnification for longest advection time available

                    # Save patches together
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5))
                        # ax.add_patch(patch)

                    # Inset image
                    axins = zoomed_inset_axes(ax, zoom, loc=4) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axins, outline=False, merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axins.contourf(XE, YE, var, cmap=cmap, levels=levels)

                    if plume: # add isohalines from surface salinity

                        # for k in xrange(salt.shape[0]):
                        #     # pdb.set_trace()
                        #     axins.contour(grid['xr'], grid['yr'], salt[k,:,:].T, [32], colors='b', linewidths=0.05, alpha=0.3)
                        axins.contour(grid.xr, grid.yr, salt.mean(axis=0).T, [32], colors='b', linewidths=1, alpha=0.8)

                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cmo.matter(norm(thisfrac))
                        thispatch.set_facecolor(color)
                        axins.add_patch(thispatch)

                    # subregion of the original image
                    axins.set_xlim(x1,x2)
                    axins.set_ylim(y1,y2)
                    plt.xticks(visible=False)
                    plt.yticks(visible=False)
                    plt.setp(axins,xticks=[],yticks=[])
                    # draw a bbox of the region of the inset axes in the parent axes and
                    # connecting lines between the bbox and the inset axes area
                    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
                    plt.draw()
                    plt.show()

                ax.set_frame_on(False)
        cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
        if log:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            fig.savefig('figures/coastconn/likelihood/interannual-log-' + season + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')
        else:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
            # if whichH=='Hall':
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            # elif whichH=='H':
            #     cb.set_label('Connection to coastal boxes in ' + str(days[j]) + ' days')
            if zoomed:
                fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + 'zoomed.png', bbox_inches='tight')
            else:
                fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')
        plt.close()
        ####


def plot_seasonal():
    '''
    Plot seasonal comparison of likelihood, either overall or just certain parts.
    '''

    cmap = cmo.speed  # 'YlGn'
    log = False
    vulnerability = True  # Whether or not to plot the orange vulnerability along the coast
    lanes = True  # to plot shipping lanes over plots
    # zoomed = True # True to add in a magnified region, for the 30 days advection timing
    whichboxes = 'galveston' # 'all', 'porta', 'galveston', 'both-zoom'
    # which boxes along coast to use for vulnerability. 0:342 for all
    # Port A:
    if whichboxes=='all':
        boxes = np.arange(0,342)
        whichH = 'Hall' # Hall for histogram for entire coastline at once; H for by coast box
        zoomed = False
    elif whichboxes=='all2':
        boxes = np.arange(0,342)
        whichH = 'H' # use coast box histograms instead of combined
        zoomed = False
    elif whichboxes=='porta':
        boxes = np.arange(103,123)
        whichH = 'H'
        zoomed = True
        # limits for zoomed box
        # if season == 'winter':
        x1, x2, y1, y2 = 86000, 340800, 465000, 715000
        #     x1, x2, y1, y2 = 86000, 277100, 527500, 715000
        # elif season == 'summer':
        #     x1, x2, y1, y2 = 86000, 277100, 465000, 652500
        zoom = 2.0
    elif whichboxes=='galveston':
        boxes = np.arange(160,180)
        whichH = 'H'
        zoomed = True
        x1, x2, y1, y2 = 277100, 531900, 660000, 810000
        zoom = 2.75
    elif whichboxes=='both-zoom':
        boxes = np.arange(103,123)  # analysis is for Port A but zoom box shows both areas
        whichH = 'H'
        zoomed = True
        x1, x2, y1, y2 = 86000, 531900, 465000, 810000
        zoom = 1.25

    d = np.load('calcs/coastpaths.npz')
    pathsxy = d['pathsxy']
    d.close()

    ## Read in files ##
    filename = 'calcs/coastconn/likelihood/hist-seasonal-' + whichboxes + '.npz'
    if not os.path.exists(filename):
        Files = []
        Files.append(glob.glob('calcs/coastconn/likelihood/hist-20??-0[1,2].npz'))
        Files.append(glob.glob('calcs/coastconn/likelihood/hist-20??-0[7,8].npz'))
        H = np.zeros((len(Files),6,100,100))
        ndbox = np.zeros((len(Files),6,342))
        for i,files in enumerate(Files): # winter and summer
            numfiles = 0
            for File in files: # months/years within winter or summer
                print File
                d = np.load(File)
                numfiles += d['numfiles']
                # size [simulations in file x advection days x coast box]
                idstemp = d['ids'].item()
                for k in xrange(d['numfiles']): # loop over simulations in file which don't want to mix up
                    for j in xrange(6): #advection days
                        # add together the drifters for the desired boxes
                        ids = []
                        [ids.extend(idstemp[k,j,box]) for box in boxes]
                        ids = list(set(ids)) # eliminate nonunique drifters
                        Htemp, _, _ = np.histogram2d(xp0[ids], yp0[ids], bins=bins,
                                        range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                        H[i,j] += Htemp
                days = d['days']
                xe = d['xe']; ye = d['ye']
                ndbox[i,:,:] += d['ndbox']
                d.close()
            # Divide by number of starting drifters
            # if (whichH=='Hall'):
            H[i] /= (numfiles*Hstart)
        np.savez(filename, H=H, days=days, xe=xe, ye=ye, ndbox=ndbox)
    else:
        d = np.load(filename)
        H = d['H']; days = d['days']; xe = d['xe']; ye = d['ye']; ndbox = d['ndbox']
        d.close()
    ####
    # pdb.set_trace()
    # limit boxes
    pathsxy = pathsxy[boxes]
    ndbox = ndbox[:,:,boxes]

    XE, YE = np.meshgrid(tracpy.op.resize(xe, 0), tracpy.op.resize(ye, 0))

    # Loop through advection days
    for j in xrange(H.shape[1]):

        ## Plot setup ##
        fig, axarr = plt.subplots(1,2)#, sharex=True)
        fig.set_size_inches(13, 6.6125)
        fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
        for i, ax in enumerate(axarr): # loop through seasons
           # Titles for subplots
            if i==0:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                ax.set_title('Winter')
            elif i==1:
                tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
                ax.set_title('Summer')
            # if whichH=='Hall':
            levels = np.linspace(0,100,11)
            var = H[i,j].T*100.
            # elif whichH=='H':
                # levels = np.linspace(0,1,11)
                # var = H[i,j].T/H[:,j].max()
            if log:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels, norm=colors.LogNorm())
            else:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels)

                if vulnerability:
                    # Add on vulnerability of coastline
                    # Need to plot the coast boxes as patches and color them according to vulnerability level
                    # http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
                    # we need to normalize the data to 0..1 for the full
                    # range of the colormap
                    fracs = ndbox[i,j,:].astype(float)/ndbox[:,j,:].max() # max across seasons
                    norm = colors.Normalize(fracs.min(), fracs.max())

                    # Save patches together
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))

                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cmo.matter(norm(thisfrac))
                        thispatch.set_facecolor(color)
                        # pdb.set_trace()
                        ax.add_patch(thispatch)

                if lanes:
                    # plot shipping lanes
                    for entry in proj.fairway:
                        entry = np.asarray(entry)
                        ax.plot(entry[:, 0], entry[:, 1], 'deepskyblue', lw=0.7, alpha=0.6)

                if zoomed: # and j==H.shape[1]-1: # magnification for longest advection time available

                    if vulnerability:
                        # Save patches together
                        patches = []
                        for path in pathsxy:
                            patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))
                            # ax.add_patch(patch)

                    # Inset image
                    axins = zoomed_inset_axes(ax, zoom, loc=4) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axins, outline=[0, 0, 0, 0], merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axins.contourf(XE, YE, var, cmap=cmap, levels=levels)

                    if vulnerability:
                        # assign shades of colormap to the patches according to values, and plot
                        for thisfrac, thispatch in zip(fracs, patches):
                            color = cmo.matter(norm(thisfrac))
                            thispatch.set_facecolor(color)
                            axins.add_patch(thispatch)

                    # subregion of the original image
                    # x1, x2, y1, y2 = 86000, 340800, 465000, 715000
                    axins.set_xlim(x1,x2)
                    axins.set_ylim(y1,y2)
                    plt.xticks(visible=False)
                    plt.yticks(visible=False)
                    plt.setp(axins,xticks=[],yticks=[])
                    # draw a bbox of the region of the inset axes in the parent axes and
                    # connecting lines between the bbox and the inset axes area
                    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
                    if lanes:
                        # plot shipping lanes
                        for entry in proj.fairway:
                            entry = np.asarray(entry)
                            axins.plot(entry[:, 0], entry[:, 1], 'deepskyblue', lw=0.7, alpha=0.6)
                    plt.draw()
                    plt.show()

            ax.set_frame_on(False)
        cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
        if log:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            fig.savefig('figures/coastconn/likelihood/seasonal-log-' + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')
        else:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
            # if whichH=='Hall':
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            # elif whichH=='H':
            #     cb.set_label('Connection to coastal boxes in ' + str(days[j]) + ' days')
            if zoomed:
                if lanes:
                    fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + 'zoomed-lanes.png', bbox_inches='tight')
                else:
                    fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + 'zoomed.png', bbox_inches='tight')
            else:
                fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')
        plt.close()
        ####


def likelihood():
    '''
    Aggregate likelihood of connection from locations with coast
    in different time periods.
    '''

    nd = np.load('calcs/xyg0.npz')['xg0'].size # # of drifters

    # Loop through along-coast boxes to find which other boxes they are connected to
    years = np.arange(2014,2015)
    months = [1,2,7,8]
    days = np.array([3,5,10,15,20,30])
    for year in years:
        for month in months:
            fname = 'calcs/coastconn/likelihood/hist-' + str(year) + '-' + str(month).zfill(2) + '.npz'
            if not os.path.exists(fname):
                Files = glob.glob('calcs/alongcoastconn/' + str(year) \
                            + '-' + str(month).zfill(2) + '-*T0*.npz')

                # likelihood histogram for each advection time examined
                Hall = np.zeros((days.size, bins[0], bins[1]))
                # number of drifters reaching each coast box (yes or no, no weighting) in day days
                ndbox = np.zeros((days.size, len(pts)))
                # http://stackoverflow.com/questions/4064277/2d-array-of-lists-in-python
                # will be size [days in month/year x advection days x coast box x drifters added for files]
                ids = defaultdict(lambda  : defaultdict(list)) # drifter ids, will have to convert to starting location later

                for k, File in enumerate(Files):
                    print File
                    d = np.load(File)
                    # [# of box paths x # drifters that enter a box x 5 (max # of crosses checked for)]
                    inbox = d['inbox'] # time in decimal days when a drifter enters a box path
                    # outbox = d['outbox'] # time in decimal days when a drifter exists a box path
                    # inds are referenced to the drifters in the shelf transport runs
                    inds = d['iinside'] # indices from the original drifters corresponding to in/outbox
                    d.close()

                    # code to switch between sets of indices
                    # this has Trues for the drifters for this simulation that enter
                    # the outer path
                    code = np.zeros(nd); code[inds] = 1; code = find(code.astype(bool))
                    # xp0[code] gives the x positions in projected space of drifters that were
                    # analyzed in the code due to entering the coastal boxes at some point.

                    # xp0temp = xp0[code]
                    # xp0temp = xp0temp[np.newaxis,:].repeat(inbox.shape[0], axis=0)
                    # yp0temp = yp0[code]
                    # yp0temp = yp0temp[np.newaxis,:].repeat(inbox.shape[0], axis=0)

                    for i, day in enumerate(days): # loop through number of advection days

                        # Drifters that enter a coast box within day days [coast box x set of drifters]
                        ind = (inbox[:,:,0]<=day)

                        # How many drifters enter each box by day days?
                        ndbox[i,:] += ind.sum(axis=1)

                        # overall histogram
                        # which drifters enter a coast box in this time period
                        tocoast = (ind.sum(axis=0)).astype(bool)
                        # xp = xp0[code][ind]; yp = yp0[code][ind]
                        # # Original xp, yp locations of drifters that enter a coast box within day days
                        # xpsaveall = xp0temp[ind]; ypsaveall = yp0temp[ind]
                        # xpsaveall = xp0[code][tocoast]; ypsaveall = yp0[code][tocoast];
                        # # add in starting locations for drifters that don't exit a coast box
                        pttemp = []
                        [pttemp.extend(pts[j]) for j in xrange(len(pts))]
                        inds2use = set(np.concatenate((code[tocoast], pttemp)))
                        Halltemp, _, _ = np.histogram2d(xp0[list(inds2use)], yp0[list(inds2use)], bins=bins,
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                        Hall[i,:,:] += Halltemp
                        # pdb.set_trace()

                        # loop through each coast box path to calculate separate origin histogram
                        for j in xrange(len(pts)):

                            # indices of drifters that start in this box, referenced to shelf transport seeding
                            pt = pts[j]

                            # Keep every simulation set of drifters separate
                            ids[k,i,j] = code[ind[j,:]] # drifter indices/ids as referenced to original seed locations

                            # add in drifters that start in box in case they didn't exit and so haven't been counted
                            ids[k,i,j] = np.concatenate((ids[k,i,j], pt))

                        # # This is how to reconstitute the histogram for some number of boxes: combine the drifter
                        # # ids together for all of the boxes, use set to find unique values, make histogram of the
                        # # unique starting locations. When used with all of the boxes in a test case, this recreated
                        # # Halltemp above.
                        # # loop over one simulation
                        # idstest = [] #np.empty(0)
                        # [idstest.extend(ids[0,0,j]) for j in xrange(342)]
                        # idstest = set(idstest) # eliminate extras
                        # # loop over another simulation
                        # idstest2 = [] #np.empty(0)
                        # [idstest2.extend(ids[1,0,j]) for j in xrange(342)]
                        # idstest2 = set(idstest2) # eliminate extras
                        # Htemp, _, _ = np.histogram2d(xp0[list(idstest)+list(idstest2)], yp0[list(idstest)+list(idstest2)], bins=bins,
                        #                 range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])

                # Save the month/year's worth of histograms
                # numfiles is to calculate the number drifters from bins for the the number of runs
                # aggregated together, compared with the appropriate number of starting drifters overall
                np.savez(fname, ids=dict(ids), xe=xe, ye=ye, days=days, numfiles=len(Files), ndbox=ndbox, Hall=Hall)


def calc_metrics():
    '''
    Calculate metrics related to oil coastal connectivity to compare with TXLA statistics
    for correlations.
    '''

    ## WINTER ##

    years = np.arange(2004, 2015)

    like = [] # integrated likelihood over domain of hitting the coast
    V = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast)
    v = [] # location of peak vulnerability
    Vporta = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast) Port A
    Vgalveston = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast) Galveston

    likefname = 'calcs/coastconn/likelihood/hist-winterinterannual.npz'
    d = np.load(likefname)
    H = d['H'] # interannual histograms, [years x 6 x 100 x 100]
    ndbox = d['ndbox'] # along-coast vulnerability [years x 6 x 342]
    d.close()

    # Loop through 2004-2014
    for i, year in enumerate(years):
        like.append(np.nansum(H[i,-1,:,:])) # use 30 day advection
        V[i,:] = ndbox[i,:,:].sum(axis=-1)
        Vporta[i,:] = ndbox[i,:,103:123].sum(axis=-1)
        Vgalveston[i,:] = ndbox[i,:,160:180].sum(axis=-1)
        # V.append(ndbox[i,-1,:].sum()) # use 30 day advection
        v.append(ndbox[i,-1,:].argmax()) # location of max

    np.savez('calcs/coastconn/likelihood/hist-wintermetrics.npz', like=like, V=V, v=v, Vporta=Vporta, Vgalveston=Vgalveston)


    ## SUMMER ##

    like = [] # integrated likelihood over domain of hitting the coast
    V = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast)
    v = [] # location of peak vulnerability
    Vporta = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast) Port A
    Vgalveston = np.zeros((years.size, 6)) # magnitude of vulnerability (integrated over coast) Galveston

    likefname = 'calcs/coastconn/likelihood/hist-summerinterannual.npz'
    d = np.load(likefname)
    H = d['H'] # interannual histograms, [years x 6 x 100 x 100]
    ndbox = d['ndbox'] # along-coast vulnerability [years x 6 x 342]
    d.close()

    # Loop through 2004-2014
    for i, year in enumerate(years):
        like.append(np.nansum(H[i,-1,:,:])) # use 30 day advection
        V[i,:] = ndbox[i,:,:].sum(axis=-1)
        Vporta[i,:] = ndbox[i,:,103:123].sum(axis=-1)
        Vgalveston[i,:] = ndbox[i,:,160:180].sum(axis=-1)
        # V.append(ndbox[i,-1,:].sum()) # use 30 day advection
        v.append(ndbox[i,-1,:].argmax()) # location of max

    np.savez('calcs/coastconn/likelihood/hist-summermetrics.npz', like=like, V=V, v=v, Vporta=Vporta, Vgalveston=Vgalveston)


if __name__ == "__main__":
    likelihood()
