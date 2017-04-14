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
# import geopandas as gpd
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


# Get and read in shipping lanes data
# get data
# info: https://www.data.boem.gov/homepg/pubinfo/repcat/arcinfo/zipped/gomr_fairways.htm
# url = 'http://www.data.boem.gov/homepg/pubinfo/repcat/arcinfo/zipped/fairway.zip'
# dirname = 'fairway'
# if not os.path.exists(dirname):
#     os.mkdir('fairway')
# if not os.path.exists('fairway.zip'):
#     os.system('wget ' + url)
# os.system('unzip -d fairway fairway.zip')

# read in shipping lanes data
proj.readshapefile('fairway/fairway', 'fairway')


# # Read in coast sensitivity data
# # Texas, under Environmental Sensitivity Index Shoreline
# # http://www.glo.texas.gov/land/land-management/gis/
# # Louisiana, under 2003>Shapefiles/ArcView 3.x project
# # http://response.restoration.noaa.gov/maps-and-spatial-data/download-esi-maps-and-gis-data.html#Louisiana
# sens = gpd.read_file('/Users/kthyng/Documents/research/papers/coastal-impacts-MPB/revisions/figures/ESI/ESI.shp')
# # bool for if extra sensitive to oil spill
# # http://response.restoration.noaa.gov/maps-and-spatial-data/example-shoreline-ranking.html
# sens10 = ['10' in sens.ESI_1[i] for i in range(len(sens))]
# # # louisiana
# # la = gpd.read_file('/Users/kthyng/Documents/research/papers/coastal-impacts-MPB/revisions/figures/Louisiana_2003_Shapefiles/AVPROJ/SHAPE/esip.shp')
# # la10 = ['10' in la.ESI[i] for i in range(len(la))]

def plot_interannual():
    '''
    Plot interannual comparison of likelihood, either overall or just certain parts.
    '''

    cmap = 'YlGn'
    log = False
    season = 'winter' # 'winter' or 'summer'
    # zoomed = True
    whichboxes = 'porta' # 'all', 'porta', 'galveston'
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
        plume = True
    elif whichboxes=='galveston':
        boxes = np.arange(160,180)
        whichH = 'H'
        zoomed = True
        x1, x2, y1, y2 = 277100, 531900, 660000, 810000
        zoom = 2.75
        # x1, x2, y1, y2 = 277100, 531900, 560000, 810000
        # zoom = 2.0

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
                print(File)
                d = np.load(File)
                Htemp = d[whichH] # 6x100x100, overall histogram
                numfiles += d['numfiles']
                if whichH=='Hall':
                    H[i] += Htemp
                elif whichH=='H':
                    H[i] = Htemp[:,boxes,:,:].sum(axis=1)
                days = d['days']
                xe = d['xe']; ye = d['ye']
                ndbox[i,:,:] += d['ndbox']
                d.close()
            # Divide by number of starting drifters
            if whichH=='Hall':
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
                pars=np.arange(20, 36, 2), outline=[0,0,0,0], parslabels=[0, 1, 0, 0])
        elif i==11:#7:
            ax.set_axis_off()
        else:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3),
                pars=np.arange(20, 36, 2), outline=[0,0,0,0],
                merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
        if i!=11:
            if whichH=='Hall':
                levels = np.linspace(0,100,11)
                var = H[i,j].T*100.
            elif whichH=='H':
                levels = np.linspace(0,1,11)
                var = H[i,j].T/H[:,j].max()
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
                    color = cm.Oranges(norm(thisfrac))
                    thispatch.set_facecolor(color)
                    ax.add_patch(thispatch)

                if zoomed:# and j==H.shape[1]-1: # magnification for longest advection time available

                    # Save patches together
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5))
                        # ax.add_patch(patch)

                    # Inset image
                    axins = zoomed_inset_axes(ax, zoom, loc=4) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axins, outline=[0,0,0,0], merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axins.contourf(XE, YE, var, cmap=cmap, levels=levels)

                    if plume: # add isohalines from surface salinity

                        # for k in xrange(salt.shape[0]):
                        #     # pdb.set_trace()
                        #     axins.contour(grid['xr'], grid['yr'], salt[k,:,:].T, [32], colors='b', linewidths=0.05, alpha=0.3)
                        axins.contour(grid.xr, grid.yr, salt.mean(axis=0).T, [32], colors='b', linewidths=1, alpha=0.8)

                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cm.Oranges(norm(thisfrac))
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
            if whichH=='Hall':
                cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            elif whichH=='H':
                cb.set_label('Connection to coastal boxes in ' + str(days[j]) + ' days')
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
    lanes = True  # to plot shipping lanes over plots
    sensitivity = False  # overlay coastal sensitivity
    # zoomed = True # True to add in a magnified region, for the 30 days advection timing
    whichboxes = 'all' # 'all', 'porta', 'galveston'
    # which boxes along coast to use for vulnerability. 0:342 for all
    # Port A:
    # if whichboxes=='all':
    #     boxes = np.arange(0,342)
    whichH = 'Hall' # Hall for histogram for entire coastline at once; H for by coast box
    #     zoomed = False
    # elif whichboxes=='all2':
    #     boxes = np.arange(0,342)
    #     whichH = 'H' # use coast box histograms instead of combined
    #     zoomed = False
    # elif whichboxes=='porta':
    boxesPA = np.arange(103,123)
    # whichH = 'H'
    zoomed = True
        # limits for zoomed box
        # if season == 'winter':
    x1PA, x2PA, y1PA, y2PA = 86000, 340800, 465000, 700000  #705000
    # x1PA, x2PA, y1PA, y2PA = 86000, 340800, 465000, 715000
        #     x1, x2, y1, y2 = 86000, 277100, 527500, 715000
        # elif season == 'summer':
        #     x1, x2, y1, y2 = 86000, 277100, 465000, 652500
    zoomPA = 2.0
    # elif whichboxes=='galveston':
    boxesG = np.arange(160,180)
    # whichH = 'H'
    # zoomed = True
    x1G, x2G, y1G, y2G = 320000, 460000, 705000, 810000
    # x1G, x2G, y1G, y2G = 277100, 531900, 660000, 810000
    zoomG = 1.9

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
                print(File)
                d = np.load(File)
                Htemp = d[whichH] # 6x100x100, overall histogram
                numfiles += d['numfiles']
                if whichH=='Hall':
                    H[i] += Htemp
                elif whichH=='H':
                    H[i] = Htemp[:,boxes,:,:].sum(axis=1)
                days = d['days']
                xe = d['xe']; ye = d['ye']
                ndbox[i,:,:] += d['ndbox']
                d.close()
            # Divide by number of starting drifters
            if whichH=='Hall':
                H[i] /= (numfiles*Hstart)
        np.savez(filename, H=H, days=days, xe=xe, ye=ye, ndbox=ndbox)
    else:
        d = np.load(filename)
        H = d['H']; days = d['days']; xe = d['xe']; ye = d['ye']; ndbox = d['ndbox']
        d.close()
    ####

    # limit boxes
    # pathsxyPA = pathsxy[boxesPA]
    # ndboxPA = ndbox[:,:,boxesPA]
    # pathsxyG = pathsxy[boxesG]
    # ndboxG = ndbox[:,:,boxesG]

    XE, YE = np.meshgrid(tracpy.op.resize(xe, 0), tracpy.op.resize(ye, 0))

    # Loop through advection days
    for j in xrange(H.shape[1]):

        if j!=H.shape[1]-1:
            continue

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
            if whichH=='Hall':
                levels = np.linspace(0,100,11)
                var = H[i,j].T*100.
            elif whichH=='H':
                levels = np.linspace(0,1,11)
                var = H[i,j].T/H[:,j].max()
            if log:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels, norm=colors.LogNorm())
            else:
                mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels)

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
                    ax.add_patch(thispatch)

                if zoomed: # and j==H.shape[1]-1: # magnification for longest advection time available

                    # Save patches together
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))
                        # ax.add_patch(patch)
                    ## Port Aransas ##
                    # Inset image
                    axinsPA = zoomed_inset_axes(ax, zoomPA, loc=4) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axinsPA, outline=[0,0,0,0], merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axinsPA.contourf(XE, YE, var, cmap=cmap, levels=levels)

                    if sensitivity:
                        sens[sens10].to_crs(proj.srs).plot(ax=axinsPA, color='b', linewidth=3)

                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cmo.matter(norm(thisfrac))
                        thispatch.set_facecolor(color)
                        axinsPA.add_patch(thispatch)

                    # subregion of the original image
                    # x1, x2, y1, y2 = 86000, 340800, 465000, 715000
                    axinsPA.set_xlim(x1PA,x2PA)
                    axinsPA.set_ylim(y1PA,y2PA)
                    plt.xticks(visible=False)
                    plt.yticks(visible=False)
                    plt.setp(axinsPA,xticks=[],yticks=[])
                    # draw a bbox of the region of the inset axes in the parent axes and
                    # connecting lines between the bbox and the inset axes area
                    mark_inset(ax, axinsPA, loc1=1, loc2=3, fc="none", ec="0.5")
                    # pdb.set_trace()

                    if lanes:
                        # plot shipping lanes
                        for entry in proj.fairway:
                            entry = np.asarray(entry)
                            axinsPA.plot(entry[:, 0], entry[:, 1], 'deepskyblue', lw=0.7, alpha=0.6)
                    ####

                    ## Galveston Bay ##
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))
                    # Inset image
                    axinsG = zoomed_inset_axes(ax, zoomG, loc=2) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axinsG, outline=[0,0,0,0], merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axinsG.contourf(XE, YE, var, cmap=cmap, levels=levels)

                    if sensitivity:
                        sens[sens10].to_crs(proj.srs).plot(ax=axinsG, color='b', linewidth=3)

                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cmo.matter(norm(thisfrac))
                        thispatch.set_facecolor(color)
                        axinsG.add_patch(thispatch)

                    # subregion of the original image
                    # x1, x2, y1, y2 = 86000, 340800, 465000, 715000
                    axinsG.set_xlim(x1G,x2G)
                    axinsG.set_ylim(y1G,y2G)
                    plt.xticks(visible=False)
                    plt.yticks(visible=False)
                    plt.setp(axinsG,xticks=[],yticks=[])
                    # draw a bbox of the region of the inset axes in the parent axes and
                    # connecting lines between the bbox and the inset axes area
                    mark_inset(ax, axinsG, loc1=1, loc2=3, fc="none", ec="0.5")

                    if lanes:
                        # plot shipping lanes
                        for entry in proj.fairway:
                            entry = np.asarray(entry)
                            axinsG.plot(entry[:, 0], entry[:, 1], 'deepskyblue', lw=0.7, alpha=0.6)
                    ####


                    plt.draw()
                    plt.show()
                    # pdb.set_trace()

                if lanes:
                    # plot shipping lanes
                    for entry in proj.fairway:
                        entry = np.asarray(entry)
                        ax.plot(entry[:, 0], entry[:, 1], 'deepskyblue', lw=0.7, alpha=0.6)

                if sensitivity:
                    sens[sens10].to_crs(proj.srs).plot(ax=ax, color='b', lw=2)

            ax.set_frame_on(False)
        cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
        if log:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            fig.savefig('figures/coastconn/likelihood/seasonal-log-' + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')
        else:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
            if whichH=='Hall':
                cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            elif whichH=='H':
                cb.set_label('Connection to coastal boxes in ' + str(days[j]) + ' days')
            if zoomed:
                if lanes:
                    if sensitivity:
                       fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + 'zoomed-2boxes-lanes-sensitivity.png', bbox_inches='tight')
                    else:
                       fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + 'zoomed-2boxes-lanes.png', bbox_inches='tight')
                else:
                    fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + 'zoomed-2boxes-lanes.png', bbox_inches='tight')
            else:
                fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight')

        plt.close()
        ###


def likelihood():
    '''
    Aggregate likelihood of connection from locations with coast
    in different time periods.
    '''

    nd = np.load('calcs/xyg0.npz')['xg0'].size # # of drifters

    # Loop through along-coast boxes to find which other boxes they are connected to
    years = np.arange(2013,2015)
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
                # Histograms of where drifters originate for advection times by coast box
                H = np.zeros((days.size, len(pts), bins[0], bins[1]))

                for File in Files:
                    print(File)
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
                        xpsaveall = xp0[code][tocoast]; ypsaveall = yp0[code][tocoast];
                        Halltemp, _, _ = np.histogram2d(xpsaveall, ypsaveall, bins=bins,
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                        Hall[i,:,:] += Halltemp
                        # pdb.set_trace()

                        # loop through each coast box path to calculate separate origin histogram
                        for j in xrange(len(pts)):

                            # indices of drifters that start in this box, referenced to shelf transport seeding
                            pt = pts[j]

                            # projected drifter origin locations that end up in this coast box in day days
                            xpsave = xp0[code][ind[j,:]]; ypsave = yp0[code][ind[j,:]]

                            # pdb.set_trace()

                            # Add in drifters that start in coast boxes
                            xpcoast = set(xp0[pt]) - set(xpsave)
                            xpsave = np.concatenate((xpsave, list(xpcoast)))
                            ypcoast = set(yp0[pt]) - set(ypsave)
                            ypsave = np.concatenate((ypsave, list(ypcoast)))

                            # Find histogram of xpsave, ypsave points for this simulation/box origin points
                            Htemp, _, _ = np.histogram2d(xpsave, ypsave, bins=bins,
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                            # aggregate number of drifters starting in different histogram bins
                            # that reach coastline for each month/year combination
                            H[i,j,:,:] += Htemp



                # Save the month/year's worth of histograms
                # numfiles is to calculate the number drifters from bins for the the number of runs
                # aggregated together, compared with the appropriate number of starting drifters overall
                np.savez(fname, H=H, xe=xe, ye=ye, days=days, numfiles=len(Files), ndbox=ndbox, Hall=Hall)
                # pdb.set_trace()


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


# if __name__ == "__main__":
#     likelihood()
