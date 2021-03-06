'''
Run plots for drifters from these simulations.
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
# import init
from datetime import datetime, timedelta
from glob import glob
import op
from matplotlib.mlab import find
from matplotlib import ticker
from matplotlib import cbook
import calendar
# from matplotlib.path import Path

# mpl.rcParams['text.usetex'] = True
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


def init(whichtime, whichtype, whichdir):
    '''
    Initialize necessary information for which type of plot/analysis we 
    are doing. 

    Currently hardwired to do shelf connectivity plots.

    Inputs:

    whichtime       Which plot time we are doing: 'weatherband', 'seasonal', 'interannual'
                    Controls the number and layout of subplots. Also controls which Files
                    are loaded in.
    whichtype       Which type of plot: 'cross' (shelf) or 'coast' (-al connectivity)
    '''

    if whichtype == 'cross':
        cmap = 'YlOrRd'
        base = 'calcs/shelfconn/'
    elif whichtype == 'coastCH':
        cmap = 'YlGn'
        base = 'calcs/coastconn/' + whichdir + '/CH/'
    elif whichtype == 'coastMX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/' + whichdir + '/MX/'
    elif whichtype == 'coastSTX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/' + whichdir + '/STX/'
    elif whichtype == 'coastNTX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/' + whichdir + '/NTX/'
    elif whichtype == 'coastLA':
        cmap = 'YlGn'
        base = 'calcs/coastconn/' + whichdir + '/LA/'
    elif whichtype == 'D2':
        cmap = 'BuPu'
        base = 'tracks/'

    #pdb.set_trace()
    Files = []
    if whichtime=='seasonal':
        # Seasonal returns Files that has two entries of lists
        Files.append(glob(base + '*-0[1-2]-*.*'))
        #Files.append(glob(base + '*-02-01*.npz'))
        Files.append(glob(base + '*-0[7-8]-*.*'))
        #Files.append(glob(base + '*-08-01*.npz'))
    elif whichtime == 'weatherband1':
        Files.append(glob(base + '2007-01-[16-19]-*.*'))
        Files.append(glob(base + '2007-01-20T04.*'))
    elif whichtime == 'interannual-winter':
        # interannual returns Files that has 12 entries of lists
        Files.append(glob(base + '2004-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2005-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2006-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2007-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2008-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2009-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2010-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2011-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2012-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2013-0[1-2]-*.*'))#npz'))
        Files.append(glob(base + '2014-0[1-2]-*.*'))#npz'))
    elif whichtime == 'interannual-summer':
        # interannual returns Files that has 12 entries of lists
        Files.append(glob(base + '2004-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2005-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2006-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2007-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2008-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2009-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2010-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2011-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2012-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2013-0[7-8]-*.*'))#npz'))
        Files.append(glob(base + '2014-0[7-8]-*.*'))#npz'))
    elif ('interannual' in whichtime) and not ('mean' in whichtime):
        month = whichtime.split('-')[-1]
        Files.append(glob(base + '2004-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2005-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2006-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2007-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2008-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2009-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2010-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2011-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2012-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2013-' + month + '-*.*'))#npz'))
        Files.append(glob(base + '2014-' + month + '-*.*'))#npz'))
    elif whichtime == 'interannual-mean':
        month = whichtime.split('-')[-1]
        Files.append(glob(base + '2004-*.*'))#npz'))
        Files.append(glob(base + '2005-*.*'))#npz'))
        Files.append(glob(base + '2006-*.*'))#npz'))
        Files.append(glob(base + '2007-*.*'))#npz'))
        Files.append(glob(base + '2008-*.*'))#npz'))
        Files.append(glob(base + '2009-*.*'))#npz'))
        Files.append(glob(base + '2010-*.*'))#npz'))
        Files.append(glob(base + '2011-*.*'))#npz'))
        Files.append(glob(base + '2012-*.*'))#npz'))
        Files.append(glob(base + '2013-*.*'))#npz'))
        Files.append(glob(base + '2014-*.*'))#npz'))
    elif ('monthly' in whichtime) and not ('mean' in whichtime):
        year = whichtime.split('-')[-1]
        Files.append(glob(base + year + '-01-*.*'))#npz'))
        Files.append(glob(base + year + '-02-*.*'))#npz'))
        Files.append(glob(base + year + '-03-*.*'))#npz'))
        Files.append(glob(base + year + '-04-*.*'))#npz'))
        Files.append(glob(base + year + '-05-*.*'))#npz'))
        Files.append(glob(base + year + '-06-*.*'))#npz'))
        Files.append(glob(base + year + '-07-*.*'))#npz'))
        Files.append(glob(base + year + '-08-*.*'))#npz'))
        Files.append(glob(base + year + '-09-*.*'))#npz'))
        Files.append(glob(base + year + '-10-*.*'))#npz'))
        Files.append(glob(base + year + '-11-*.*'))#npz'))
        Files.append(glob(base + year + '-12-*.*'))#npz'))
    elif whichtime == 'monthly-mean':
        Files.append(glob(base + '20??-01-*.*'))#npz'))
        Files.append(glob(base + '20??-02-*.*'))#npz'))
        Files.append(glob(base + '20??-03-*.*'))#npz'))
        Files.append(glob(base + '20??-04-*.*'))#npz'))
        Files.append(glob(base + '20??-05-*.*'))#npz'))
        Files.append(glob(base + '20??-06-*.*'))#npz'))
        Files.append(glob(base + '20??-07-*.*'))#npz'))
        Files.append(glob(base + '20??-08-*.*'))#npz'))
        Files.append(glob(base + '20??-09-*.*'))#npz'))
        Files.append(glob(base + '20??-10-*.*'))#npz'))
        Files.append(glob(base + '20??-11-*.*'))#npz'))
        Files.append(glob(base + '20??-12-*.*'))#npz'))

    # pdb.set_trace()

    return Files, cmap


def calc_histogram(xp, yp, whichtype, bins=(60,60), 
                    Xrange=None, Yrange=None):
    '''
    Calculate the histograms for the connectivity calculations.
    Calculate both the histogram of starting locations and of connected 
    trajectories.

    Inputs:

    xp, yp          Drifter x,y positions in projected space
    bins            ([60,60]) x and y numbers of bins for the histogram
    Xrange          (None) min and max values for x direction in format (xmin, xmax)
    Yrange          (None) min and max values for y direction in format (ymin, ymax)
    '''

    if whichtype == 'cross' or 'coast' in whichtype:
        # Then we are finding the number of drifters starting in each
        # bin for division into a probability
        # pdb.set_trace()
        if Xrange!=None and Yrange!=None:
            H, xe, ye = np.histogram2d(xp.flatten(), yp.flatten(), bins=bins, 
                                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
        else:
            H, xe, ye = np.histogram2d(xp, yp, bins=bins)

    elif whichtype == 'D2' or whichtype == 'fsle':
        # We are finding which drifters started in each bin, by index
        xes = np.linspace(Xrange[0], Xrange[1], bins[0])
        yes = np.linspace(Yrange[0], Yrange[1], bins[1])
        # Initialize an array of lists: http://mail.scipy.org/pipermail/numpy-discussion/2009-November/046566.html
        filler = np.frompyfunc(lambda x: list(), 1, 1)
        H = np.empty((yes.size-1,xes.size-1), dtype=np.object)
        filler(H, H)
        for i, xe in enumerate(xes[:-1]): # loop through edges in x
            for j, ye in enumerate(yes[:-1]): # loop through edges in y
                # H contains indices of corresponding drifter seed locations
                # THESE WILL NEED TO BE LISTS IN THE ARRAY
                # pdb.set_trace()
                H[j,i].append(find((xp>xe) * (xp<xes[i+1]) * (yp>ye) * (yp<yes[j+1])))
        # pdb.set_trace()
        xe = xes; ye = yes;

    return H, xe, ye


def calc_metric(xp, yp, Hstart, whichtype, r=[0,1.05]):
    '''
    Calculate metric given by whichtype 
    '''

    # loop through histogram bins and calculate metrics for each bin
    # pdb.set_trace()
    metric = np.empty((Hstart.shape[0], Hstart.shape[1], xp.shape[1]))
    nnans = np.empty((Hstart.shape[0], Hstart.shape[1], xp.shape[1]))
    for j in xrange(Hstart.shape[0]):
        for i in xrange(Hstart.shape[1]):

            if whichtype == 'D2':

                # D^2 is already averaged
                # Picks out the drifters (xp,yp) that are in the j,i-th histogram bin, and looks at their statistics
                # metric[j,i,:], nnans[j,i,:], pairs = tracpy.calcs.rel_dispersion(xp[Hstart[j,i][0],:], yp[Hstart[j,i][0],:], r=1.05, squared=False)
                metric[j,i,:], nnans[j,i,:], pairs = tracpy.calcs.rel_dispersion(xp[Hstart[j,i][0],:], yp[Hstart[j,i][0],:], r=r, squared=True)
                # np.savez('calcs/pairs/bin' + str(i) + '_' + str(H.size) + '.npz')
            
            elif whichtype == 'fsle':
                tSavetemp = tracpy.calcs.calc_fsle(lonp, latp, tp, alpha=np.sqrt(2))
                ind = ~np.isnan(tSavetemp)
                metric = np.nansum(tSavetemp, axis=0)
                nnans = ind.sum(axis=0)
                # still have to calculate fsle
                # metric = 1./(tSave/nnans) # CHECK THIS

    return metric, nnans


def plot_setup(whichtime, grid):
    '''
    Set up the figure for plotting based on the string in which.

    Inputs:

    whichtime           Which type of plot to do
    grid            Previously-read-in grid dictionary, with basemap included
    '''

    if whichtime=='seasonal':
        # fig = plt.figure(figsize=(22,10))
        fig, axarr = plt.subplots(1,2)#, sharex=True)
        # pdb.set_trace()
        fig.set_size_inches(13, 6.6125)
        fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
        # ax = fig.add_subplot(1,2,1)
        for i, ax in enumerate(axarr):
           # Titles for subplots
            if i==0:
                 tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), col='0.3', halpha=0.5)
                 ax.set_title('Winter')
                 # ax.text(0.45, 0.95, 'Winter', transform=ax.transAxes)#, bbox=dict(facecolor='white', edgecolor=None), fontsize=14)#, alpha=0.5))
            elif i==1:
                 tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2), col='0.3', halpha=0.5)
                 ax.set_title('Summer')
            # suptitle
            #fig.suptitle('Probability of material crossing the shelf in 30 days, 2004-2010', y=.94)
            ax.set_frame_on(False)

    elif 'interannual' in whichtime: # summer or winter

        fig, axarr = plt.subplots(4,3)
        fig.set_size_inches(8.9, 11.5)
        fig.subplots_adjust(left=0.04, bottom=0.1, right=1.0, top=0.99, wspace=0.0001, hspace=0.05)

        # fig, axarr = plt.subplots(2,4)
        # fig.set_size_inches(13.4, 6.6125)
        # fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

        for i, ax in enumerate(axarr.flatten()):
            if i==0 or i==3 or i==6:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], parslabels=[1, 0, 0, 0],
                    merslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)
                if i==0:  # plot analysis box
                    lonbox = [-90, -91, -92, -93, -94, -95, -96, -97, -97.5, 
                              -97.5, -97, -96, -95, -94, -93, -92, -91, -90, -90]
                    latbox = [27.5, 27.5, 27.5, 27.5, 27.5, 27.5, 27.5, 27.5, 27.5, 
                              29.9, 29.9, 29.9, 29.9, 29.9, 29.9, 29.9, 29.9, 29.9, 27.5]
                    xbox, ybox = grid['basemap'](lonbox, latbox)
                    ax.plot(xbox, ybox, '0.2', lw=1)
            elif i==9:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], merslabels=[0, 0, 0, 1],
                    parslabels=[1, 0, 0, 0], col='0.3', halpha=0.5)
            elif i==10:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], merslabels=[0, 0, 0, 1],
                    parslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)
            elif i==11:
                ax.set_axis_off()
            else:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0],
                    merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)

            if i!=11:
                ax.text(0.07, 0.88, str(2004+i), transform=ax.transAxes)

            ax.set_frame_on(False)

    elif 'monthly' in whichtime:

        fig, axarr = plt.subplots(4,3)
        fig.set_size_inches(8.7, 11.5)
        fig.subplots_adjust(left=0.008, bottom=0.1, right=1.0, top=0.98, wspace=0.005, hspace=0.1)

        for i, ax in enumerate(axarr.flatten()):
            if i==0 or i==3 or i==6:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], parslabels=[1, 0, 0, 0],
                    merslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)
            elif i==9:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], merslabels=[0, 0, 0, 1],
                    parslabels=[1, 0, 0, 0], col='0.3', halpha=0.5)
            elif i==10 or i==11:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0], merslabels=[0, 0, 0, 1],
                    parslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)
            else:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 3), 
                    pars=np.arange(20, 36, 2), outline=[0,0,0,0],
                    merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0], col='0.3', halpha=0.5)

            ax.text(0.07, 0.88, calendar.month_name[i+1], transform=ax.transAxes)

            ax.set_frame_on(False)

    return fig, axarr

def plot_stuff(xe, ye, H, cmap, grid, shelf_depth, ax, levels=np.linspace(0,100,11), extend='max'):
    '''
    Do the main plotting stuff.
    '''

    XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))

    # Try with pcolor too
    #pdb.set_trace()
    mappable = ax.contourf(XE, YE, H, cmap=cmap, levels=levels, extend=extend)
    ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

    return mappable


def plot_colorbar(fig, mappable, whichtype, ticks=None, whichdir='forward', whichtime='seasonal'):
    '''
    Add colorbar to figure.
    '''

    # Make colorbar
    # Horizontal colorbar below plot
    if whichtime=='seasonal':
        cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    elif ('interannual' in whichtime) or ('monthly' in whichtime):
        cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')

    if whichtype == 'cross':
        cb.set_label('Probability of drifters crossing shelf (%)')
    elif 'coast' in whichtype: 
        cb.set_label('Probability of drifters reaching coastline region (%)')
    elif whichtype == 'D2':
        cb.set_label('Mean squared separation distance [km$^2\!$]', fontsize=13.5)
    elif whichtype == 'diff':
        cb.set_label('Probability difference [%]')

    # this trumps other settings
    if whichdir == 'back':
        cb.set_label('Normalized drifter occurrence')

    if not ticks==None:
        cb.set_ticks(ticks)


def plot_finish(fig, whichtype, whichtime, shelf_depth, itind, r, numdays):
    '''
    Save and close figure
    '''

    if whichtype == 'cross':
        fname = 'figures/' + whichtype + '/' + whichtime + str(shelf_depth) + 'advectiondays' + str(numdays) + '.png'
    elif 'coast' in whichtype: 
        fname = 'figures/' + whichtype + '/' + whichtime + '.png'
    elif whichtype == 'D2':
        fname = 'figures/' + whichtype + '/r' + str(int(r[1])) + '/' + whichtime + str(itind) + '.png'

    fig.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_diff():
    '''
    Plot difference in transport between two histograms. Plot separately for convenience,
    since this uses information previously calculated by other functions.
    import make_plots
    make_plots.plot_diff()
    '''

    # Which timing of plot: 'interannual', 'seasonal'; 'interannual-summer', 'interannual-winter'
    whichtime = 'seasonal'
    # Which type of plot: 'cross'; 'mean' (difference from interannual mean)
    whichtype = 'cross'

    shelf_depth = 100

    # Grid info
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # grid = tracpy.inout.readgrid(loc, usebasemap=True)
    # grid_filename = '../../grid.nc'
    # grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)
    grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
    grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

    if whichtype == 'cross':
        cmap = 'RdBu'
    elif whichtype == 'mean':
        cmap = 'PRGn'

    if whichtime == 'seasonal':
        # d = np.load('calcs/shelfconn/seasonal' + str(shelf_depth) + 'H.npz')
        d = np.load('figures/cross/seasonal' + str(shelf_depth) + 'H.npz')
        Hboth = d['H']; xe=d['xe']; ye=d['ye']
        H = Hboth[0,:] - Hboth[1,:]
        d.close()
    elif whichtime == 'interannual':
        dw = np.load('figures/cross/interannual-winter' + str(shelf_depth) + 'H.npz')
        Hw = dw['H']; xe=dw['xe']; ye=dw['ye']
        dw.close()
        ds = np.load('figures/cross/interannual-summer' + str(shelf_depth) + 'H.npz')
        Hs = ds['H']; xe=ds['xe']; ye=ds['ye']
        ds.close()
        H = Hw - Hs
    elif whichtime == 'interannual-summer':
        ds = np.load('figures/cross/interannual-summer' + str(shelf_depth) + 'H.npz')
        Hs = ds['H']; xe=ds['xe']; ye=ds['ye']
        ds.close()
        # pdb.set_trace()
        H = Hs - Hs.mean(axis=0)
    elif whichtime == 'interannual-winter':
        ds = np.load('figures/cross/interannual-winter' + str(shelf_depth) + 'H.npz')
        Hs = ds['H']; xe=ds['xe']; ye=ds['ye']
        ds.close()
        # pdb.set_trace()
        H = Hs - Hs.mean(axis=0)

    # if (shelf_depth==100) and (whichtime == 'seasonal'):
    #     H *= 100 


    if whichtime == 'seasonal':
        if shelf_depth == 20:
            levels = np.arange(-65,75,10)
            ticks = [-65,-45,-25,0,25,45,65]
        elif shelf_depth == 50:
            levels = np.arange(-55,65,10)
            ticks = [-55,-35,-15,0,15,35,55]
        elif shelf_depth == 100:
            levels = np.arange(-45,51,6)
            ticks = [-39,-27,-15,0,15,27,39]
            # levels = np.arange(-39,45,6)
            # ticks = [-39,-27,-15,0,15,27,39]
        fig = plt.figure(figsize=(6.8375, 6.6125))
        fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)
        ax = fig.add_subplot(111)
        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), col='0.3', halpha=0.5)
        ax.set_title('Winter-Summer Transport')
        mappable = plot_stuff(xe, ye, H.T, cmap, grid, shelf_depth, ax, levels=levels, extend='neither')
        # XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))
        # cs = ax.contour(XE, YE, H.T, [-21,21], colors='b')
        ax.set_axis_off()

        plot_colorbar(fig, mappable, 'diff', ticks=ticks, whichtime=whichtime)
        fig.text(0.125, 0.075, 'Summer', color='#cc0027')
        fig.text(0.760, 0.075, 'Winter', color='#1b72b7')


    elif 'interannual' in whichtime:

        fig, axarr = plot_setup('interannual', grid)

        if whichtype == 'mean':
            if 'summer' in whichtime:
                if shelf_depth == 20:
                    levels = np.arange(-45,48,6)
                    ticks = [-45,-33,-21,-9,0,9,21,33,45]
                elif shelf_depth == 50:
                    levels = np.arange(-52,60,8)
                    ticks = [-44,-28,-12,0,12,28,44]
                elif shelf_depth == 100:
                    levels = np.arange(-60,68,8)
                    ticks = [-60,-44,-28,-12,0,12,28,44,60]
            elif 'winter' in whichtime:
                if shelf_depth == 20:
                    levels = np.arange(-52,60,8)
                    ticks = [-44,-28,-12,0,12,28,44]
                elif shelf_depth == 50:
                    levels = np.arange(-68,76,8)
                    ticks = [-60,-44,-28,-12,0,12,28,44,60]
                elif shelf_depth == 100:
                    levels = np.arange(-60,68,8)
                    ticks = [-60,-44,-28,-12,0,12,28,44,60]

        elif whichtype == 'cross':
            if shelf_depth == 20:
                levels = np.arange(-95,105,10)
                ticks = [-85,-65,-45,-25,0,25,45,65,85]
            elif shelf_depth == 50:
                levels = np.arange(-95,105,10)
                ticks = [-85,-65,-45,-25,0,25,45,65,85]
            elif shelf_depth == 100:
                levels = np.arange(-85,95,10)
                ticks = [-85,-65,-45,-25,0,25,45,65,85]

        for i in xrange(H.shape[0]): # 1 for each subplot
            mappable = plot_stuff(xe, ye, H[i,:,:].T, cmap, grid, shelf_depth, axarr.flatten()[i], 
                                    levels=levels, extend='neither')
            
        plot_colorbar(fig, mappable, 'diff', ticks=ticks, whichtime=whichtime)

        # Labels for colorbar ends
        if ('summer' in whichtime) or ('winter' in whichtime):
            fig.text(0.125, 0.075, 'Less than mean', color='#750077')
            fig.text(0.760, 0.075, 'More than mean', color='#007432')
        else:
            fig.text(0.18, 0.075, 'Summer', color='#cc0027')
            fig.text(0.760, 0.075, 'Winter', color='#1b72b7')

    fig.savefig('figures/cross/' + whichtime + 'diff' + str(shelf_depth) + '.png', dpi=300)
    plt.close()


def run():

    # Which timing of plot: 'weatherband[1-3]', 'seasonal', 'interannual-winter', 'interannual-summer'
    # 'interannual-01' through 'interannual-12', 'monthly-2004' through 'monthly-2014'
    # 'interannual-mean' 'monthly-mean'
    whichtime = 'seasonal'
    # Which type of plot: 'cross', 'coastCH', 'coastMX', 'coastLA', 
    #  'coastNTX', 'coastSTX', 'fsle', 'D2'
    whichtype = 'cross'
    # 'forward' or 'back' in time
    whichdir = 'forward'

    # if doing D2, choose the initial separation distance too. And give a little extra for roundoff and projections
    # Also, now given as an interval.
    r = [4.75, 5.25] # [0.95, 1.05] # [4.75, 5.25] # 5.25 and 1.05

    #levels = np.linspace(0,1,11)

    shelf_depth = 100 # do 100 50 and 20 
    ishelf_depth = 2 # 2 1 0 index in cross array
    numdays = 15  # 30. Number of analysis days to consider

    # Whether to overlay previously-calculated wind stress arrows
    # from projects/txla_plots/plot_mean_wind.py on Rainier
    addwind = 0  # for adding the wind on
    years = np.arange(2004,2015) # this is just for the wind I think

    # Number of bins to use in histogram
    bins = (100,100) #(30,30)

    # Load in Files to read from based on which type of plot is being run
    Files, cmap = init(whichtime, whichtype, whichdir)

    # Grid info
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # loc = '/Users/kthyng/Documents/research/postdoc/grid.nc'
    grid = tracpy.inout.readgrid(loc, usebasemap=True)
    # grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    # vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
    # grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)
    # grid_filename = '../../grid.nc'
    # grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)

    if whichtype == 'D2':
        Hfilename = 'figures/' + whichtype + '/r' + str(int(r[1])) + '/' + whichtime + 'H.npz'
    elif 'coast' in whichtype:
        Hfilename = 'figures/' + whichtype + '/' + whichtime + 'H.npz'
    elif whichtype == 'cross':
        Hfilename = 'figures/' + whichtype + '/' + whichtime + str(shelf_depth) + 'advection' + str(numdays) + 'days-' + 'H.npz'
        # Hfilename = 'calcs/shelfconn/' + whichtime + str(shelf_depth) + 'H.npz'

    if not os.path.exists(Hfilename): 

        ## Calculate starting position histogram just once ##
        # Read in connectivity info (previously calculated). 
        # Drifters always start in the same place.
        # pdb.set_trace()
        if (whichtype == 'cross') or ('coast' in whichtype and whichdir == 'forward'):
            d = np.load(Files[0][0])

            # Calculate xrange and yrange for histograms
            Xrange = [grid['xpsi'].min(), grid['xpsi'].max()]
            Yrange = [grid['ypsi'].min(), grid['ypsi'].max()]

            # Histogram of starting locations
            if whichtype == 'cross': # results are in xg, yg
                xp0, yp0, _ = tracpy.tools.interpolate2d(d['xg0'], d['yg0'], grid, 'm_ij2xy')
            elif 'coast' in whichtype:  # results are in xp, yp
                xp0 = d['xp0']; yp0 = d['yp0']

            d.close()

        elif ('coast' in whichtype and whichdir == 'back'):

            # Calculate xrange and yrange for histograms
            Xrange = [grid['xpsi'].min(), grid['xpsi'].max()]
            Yrange = [grid['ypsi'].min(), grid['ypsi'].max()]

        elif whichtype == 'D2': # results are in xg, yg

            # Histogram of starting locations
            Hstartfile = 'calcs/dispersion/hist/Hstart_bins' + str(bins[0]) + '.npz'

            if not os.path.exists(Hstartfile): # just read in info
                d = netCDF.Dataset(Files[0][0])

                # Calculate xrange and yrange for histograms in ll
                Xrange = [grid['lonpsi'].min(), grid['lonpsi'].max()]
                Yrange = [grid['latpsi'].min(), grid['latpsi'].max()]
            
                # Histogram of starting locations
                # xp0, yp0 are lonp, latp in this case
                xp0, yp0, _ = tracpy.tools.interpolate2d(d.variables['xg'][:,0], d.variables['yg'][:,0], grid, 'm_ij2ll')

                d.close()

        if 'Hstartfile' in locals() and os.path.exists(Hstartfile): # just read in info
            Hstartf = np.load(Hstartfile)
            xe = Hstartf['xe']; ye = Hstartf['ye']
            Hstart = Hstartf['Hstart']
        # elif 'coast' in whichtype and whichdir == 'back':
        #     continue
        else:
            if whichdir == 'back': # aren't dividing at the end by the starting number
                Hstart = np.ones(bins)
            else:
                # For D2 and fsle, Hstart contains indices of drifters seeded in bins
                Hstart, xe, ye = calc_histogram(xp0, yp0, whichtype, bins=bins, Xrange=Xrange, Yrange=Yrange)

                if whichtype == 'D2':
                    xe, ye = grid['basemap'](xe, ye) # change from lon/lat
                    np.savez(Hstartfile, Hstart=Hstart, xe=xe, ye=ye) 

        # For D2 and fsle, H contains the metric calculation averaged over that bin
        if whichtype == 'D2' or whichtype == 'fsle':
            H = np.zeros((len(Files), Hstart.shape[0], Hstart.shape[1], 901))
            nnans = np.zeros((len(Files), Hstart.shape[0], Hstart.shape[1], 901))
        elif ('coast' in whichtype) or (whichtype=='cross'):
            H = np.zeros((len(Files), Hstart.shape[0], Hstart.shape[1]))
            nnans = np.zeros((len(Files), Hstart.shape[0], Hstart.shape[1]))

        # Loop through calculation files to calculate overall histograms
        # pdb.set_trace()
        for i, files in enumerate(Files): # Files has multiple entries, 1 for each subplot

            if whichtype == 'cross' or 'coast' in whichtype:
                Hcross = np.zeros(bins) # initialize
                #pdb.set_trace()
                HstartUse = Hstart*len(files) # multiply to account for each simulation


            for File in files: # now loop through the files for this subplot
                print File

                if whichtype == 'cross': # results are in xg, yg
                # [number of depths,number of tracks] to store time of crossing or nan if it doesn't cross
                    # Read in info
                    d = np.load(File)
                    xg0 = d['xg0']; yg0 = d['yg0']
                    cross = d['cross']
                    ind = (~np.isnan(cross[ishelf_depth,:])) * (cross[ishelf_depth, :] < numdays)
                    # ind = ~np.isnan(cross[ishelf_depth,:])  # this is for 30 days (the whole run time)
                    d.close()
                    xp, yp, _ = tracpy.tools.interpolate2d(xg0[ind], yg0[ind], grid, 'm_ij2xy')
                elif ('coast' in whichtype) and (whichdir=='forward'):  # results are in xp, yp
                    # Read in info
                    d = np.load(File)
                    xp = d['xp0']; yp = d['yp0']
                    conn = d['conn'] 
                    ind = ~np.isnan(conn)
                    xp = xp[ind]; yp = yp[ind] # pick out the drifters that reach the coastline
                    d.close()
                elif ('coast' in whichtype) and (whichdir=='back'):  # results are in xp, yp
                    # Read in info
                    d = np.load(File)
                    # backward tracks in time of drifters starting in the specified coastal area
                    xp = d['xp0']; yp = d['yp0']
                    conn = d['conn'] # indices of the drifters that started in the zone
                    # can't send nan's to histogram calculation
                    ind = ~np.isnan(xp)
                    xp = xp[ind]; yp = yp[ind] # pick out the drifters that reach the coastline
                    d.close()
                elif whichtype == 'D2' or whichtype == 'fsle':
                    sfile = 'calcs/dispersion/hist/D2/r' + str(int(r[1])) + '/' + File.split('/')[-1][:-5] + '_bins' + str(bins[0]) + '.npz'
                    if os.path.exists(sfile): # just read in info
                        already_calculated = 1
                    else:
                        already_calculated = 0
                    # # This is for distributing the workload to different processors
                    # if not '2004' in File:
                    #     continue

                    if not already_calculated:
                        print 'working on D2 ', sfile
                        d = netCDF.Dataset(File)
                        xg = d.variables['xg'][:]; yg = d.variables['yg'][:]
                        # eliminate entries equal to -1
                        ind = xg==-1
                        xg[ind] = np.nan; yg[ind] = np.nan
                        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')
                        d.close()

                # Count the drifters for the shelf_depth that have a non-nan entry
                # pdb.set_trace() 
                if whichtype == 'cross' or 'coast' in whichtype:
                    # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                    Hcrosstemp, xe, ye = calc_histogram(xp, yp, whichtype, bins=bins, Xrange=Xrange, Yrange=Yrange)
                    Hcross = np.nansum( np.vstack((Hcross[np.newaxis,:,:], Hcrosstemp[np.newaxis,:,:])), axis=0)
                elif whichtype == 'D2' or whichtype == 'fsle':
                    if not already_calculated:
                        # Calculate the metric in each bin and combine for all files
                        metric_temp, nnanstemp = calc_metric(xp, yp, Hstart, whichtype, r=r)
                        # Save calculations by bins for each file
                        # pdb.set_trace()
                        np.savez(sfile, D2=metric_temp, nnans=nnanstemp) 
                        print 'saving D2 file ', sfile
                        # metric_temp is in time, but want to show a single value for each bin in space.
                        # Take the value at the final time.
                        # pdb.set_trace()
                    else:
                        d = np.load(sfile)
                        metric_temp = d['D2']; nnanstemp = d['nnans']
                        # pdb.set_trace()
                        # filter out boxes with too few available drifters
                        ind = nnanstemp<30
                        metric_temp[ind] = np.nan

                    H[i,:] = np.nansum( np.vstack((H[np.newaxis,i,:,:,:],metric_temp[np.newaxis,:,:,:]*nnanstemp[np.newaxis,:,:,:])), axis=0) # need to un-average before combining
                    # H[i,:] = H[i,:] + metric_temp[:,:,-1]*nnanstemp[:,:,-1] # need to un-average before combining
                    nnans[i,:] = nnans[i,:] + nnanstemp[:,:,:] # need to un-average before combining

            # Calculate overall histogram
            if whichtype == 'cross' or 'coast' in whichtype:
                H[i,:] = (Hcross/HstartUse)*100
            elif whichtype == 'D2':
                # xe, ye = grid['basemap'](xe, ye) # change from lon/lat
                H[i,:] = H[i,:]/nnans[i,:]
                # np.savez('calcs/dispersion/hist/' + File.split('/')[-1][:-5] + '_bins' + str(bins[0])) 
            elif whichtype == 'fsle':
                H[i,:] = 1./H[i,:]/nnans[i,:]

        # save H
        if not os.path.exists('figures/' + whichtype): 
            os.makedirs('figures/' + whichtype)

        np.savez(Hfilename, H=H, xe=xe, ye=ye)

    else: # H has already been calculated

        Hfile = np.load(Hfilename)
        H = Hfile['H']; xe = Hfile['xe']; ye = Hfile['ye']
        #levels = np.linspace(0, np.nanmax(H), 11)

    # which time index to plot
    # a number or 'mean' or 'none' (for coast and cross)
    if whichtype == 'D2':
        itind = 100 # 100 
    elif (whichtype == 'cross') or ('coast' in whichtype):
        itind = 'none'
    # Choose consistent levels to plot
    locator = ticker.MaxNLocator(11)
    locator.create_dummy_axis()
    # don't use highest max since everything is washed out then
    # pdb.set_trace()
    # 12000 for mean interannual-summer, 20000 for mean, interannual-winter, 1400 for 100 seasonal
    # 1800 for 100 interannual-winter, 1800 for 100 interannual-summer
    if whichtype == 'D2':
        if itind == 30:
            locator.set_bounds(0, 10)
        elif itind == 100: 
            locator.set_bounds(0, 160) 
        elif itind == 150: 
            locator.set_bounds(0, 450) 
        elif itind == 300: 
            locator.set_bounds(0, 2200) 
        elif itind == 600: 
            locator.set_bounds(0, 8000) 
        elif itind == 900: 
            locator.set_bounds(0, 15000) 
        # locator.set_bounds(0, 0.2*np.nanmax(H[:,:,:,itind]))
        #locator.set_bounds(0, 0.75*np.nanmax(np.nanmax(H[:,:,:,itind], axis=1), axis=1).mean())
        levels = locator()
    elif 'coast' in whichtype and whichdir == 'back':
        hist, bin_edges = np.histogram(H.flat, bins=100) # find # of occurrences of histogram bin values
        n = np.cumsum(hist)
        Hmax = bin_edges[find(n<(n.max()-n.min())*.7+n.min())[-1]] # take the 80% of histogram occurrences as the max instead of actual max since too high
        locator.set_bounds(0, 1) 
        levels = locator()
        extend = 'max'
        H = H/Hmax
    else:
        extend = 'neither'


    # Set up overall plot, now that everything is calculated
    fig, axarr = plot_setup(whichtime, grid) # depends on which plot we're doing

    # Loop through calculation files to calculate overall histograms
    # pdb.set_trace()
    for i in xrange(H.shape[0]): # Files has multiple entries, 1 for each subplot

        # Do subplot
        # pdb.set_trace()
        # which time index to plot?
        #itind = 100
        if cbook.is_numlike(itind): # plot a particular time
            mappable = plot_stuff(xe, ye, H[i,:,:,itind], cmap, grid, shelf_depth, axarr.flatten()[i], levels=levels)
        elif itind=='mean': # plot the mean over time
            mappable = plot_stuff(xe, ye, np.nansum(H[i,:,:,:], axis=-1)/np.sum(~np.isnan(H[i,:,:,:]), axis=-1), cmap, grid, shelf_depth, axarr.flatten()[i], levels=levels)
        elif itind=='none': # just plot what is there
            if 'levels' in locals():
                mappable = plot_stuff(xe, ye, H[i,:,:].T, cmap, grid, shelf_depth, axarr.flatten()[i], extend=extend, levels=levels)
            else:
                mappable = plot_stuff(xe, ye, H[i,:,:].T, cmap, grid, shelf_depth, axarr.flatten()[i], extend=extend)
        #axarr.flatten()[i].set_title(np.nanmax(H[i,:,:,itind]))
        # Add coastline area if applicable
        if 'coast' in whichtype:
            coastloc = whichtype.split('coast')[-1]
            pts = np.load('calcs/' + coastloc + 'pts.npz')[coastloc]
            axarr.flatten()[i].plot(pts[:,0], pts[:,1], color='0.0', lw=3)
            # verts = np.vstack((pts[:,0], pts[:,1]))
            # # Form path
            # path = Path(verts.T)
            # if not path.contains_point(np.vstack((xp[jd,it],yp[jd,it]))):

        # Overlay mean wind arrows
        if addwind:
            # Right now is just for cross, interannual, winter
            year = years[i]
            # year = File.split('/')[-1].split('-')[0]
            season = whichtime.split('-')[-1]
            wind = np.load('../txla_plots/calcs/wind_stress/1st/jfm/' + str(year) + season +  '.npz')
            x = wind['x']; y = wind['y']; u = wind['u']; v = wind['v']
            q = axarr.flatten()[i].quiver(x, y, u, v, color = '0.3',
                        pivot='middle', zorder=1e35, width=0.003)
                        # scale=1.0/scale, pivot='middle', zorder=1e35, width=0.003)

            # if year == 2008:
            #     plt.quiverkey(q, 0.85, 0.07, 0.1, label=r'0.1 N m$^{2}$', coordinates='axes')



    # Add colorbar
    plot_colorbar(fig, mappable, whichtype, whichdir=whichdir, whichtime=whichtime)
    # pdb.set_trace()

    # save and close
    plot_finish(fig, whichtype, whichtime, shelf_depth, itind, r, numdays)


if __name__ == "__main__":
    run()     
