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
import init
from datetime import datetime, timedelta
from glob import glob
import op

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


def init(whichtime, whichtype):
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
        base = 'calcs/coastconn/CH/'
    elif whichtype == 'coastMX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/MX/'
    elif whichtype == 'coastSTX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/STX/'
    elif whichtype == 'coastNTX':
        cmap = 'YlGn'
        base = 'calcs/coastconn/NTX/'
    elif whichtype == 'coastLA':
        cmap = 'YlGn'
        base = 'calcs/coastconn/LA/'
    elif whichtype == 'D2':
        cmap = 'YlGnBu'
        base = 'tracks/'

    #pdb.set_trace()
    Files = []
    if whichtime=='seasonal':
        # Seasonal returns Files that has two entries of lists
        Files.append(glob(base + '*-0[1-2]-*.npz'))
        #Files.append(glob(base + '*-02-01*.npz'))
        Files.append(glob(base + '*-0[7-8]-*.npz'))
        #Files.append(glob(base + '*-08-01*.npz'))
    elif whichtime == 'weatherband1':
        Files.append(glob(base + '2007-01-[16-19]-*.npz'))
        Files.append(glob(base + '2007-01-20T04.npz'))
    elif whichtime == 'interannual-winter':
        # interannual returns Files that has 12 entries of lists
        Files.append(glob(base + '2004-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2005-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2006-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2007-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2008-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2009-0[1-2]-*.'))#npz'))
        Files.append(glob(base + '2010-0[1-2]-*.'))#npz'))
    elif whichtime == 'interannual-summer':
        # interannual returns Files that has 12 entries of lists
        Files.append(glob(base + '2004-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2005-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2006-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2007-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2008-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2009-0[7-8]-*.'))#npz'))
        Files.append(glob(base + '2010-0[7-8]-*.'))#npz'))


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
        if Xrange!=None and Yrange!=None:
            H, xe, ye = np.histogram2d(xp, yp, bins=bins, 
                                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
        else:
            H, xe, ye = np.histogram2d(xp, yp, bins=bins)

    elif whichtype == 'D2' or whichtype == 'fsle':
        # We are finding which drifters started in each bin, by index
        xes = np.linspace(Xrange[0], Xrange[1], bins[0])
        yes = np.linspace(Yrange[0], Yrange[1], bins[1])
        H = np.empty((yes.size-1,xes.size-1))
        for i, xe in enumerate(xes[:-1]): # loop through edges in x
            for j, ye in enumerate(yes[:-1]): # loop through edges in y
                # H contains indices of corresponding drifter seed locations
                # THESE WILL NEED TO BE LISTS IN THE ARRAY
                H[j,i] = xp<xe and xp>xes[i+1] and yp<ye and yp>yes[j+1]

    return H, xe, ye


def calc_metric(xp, yp, Hstart, whichtype):
    '''
    Calculate metric given by whichtype 
    '''

    if whichtype == 'D2':
        metric, nnans = tracpy.calcs.rel_dispersion(xp, yp, r=1, squared=True)
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
        fig.set_size_inches(13.675, 6.6125)
        fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)
        # ax = fig.add_subplot(1,2,1)
        for i, ax in enumerate(axarr):
           # Titles for subplots
            if i==0:
                 tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                 ax.set_title('Winter')
            elif i==1:
                 tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
                 ax.set_title('Summer')
            # suptitle
            #fig.suptitle('Probability of material crossing the shelf in 30 days, 2004-2010', y=.94)
    elif 'interannual' in whichtime: # summer or winter

        fig, axarr = plt.subplots(2,4)
        fig.set_size_inches(13.4, 6.6125)
        fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

        for i, ax in enumerate(axarr.flatten()):
           # Titles for subplots
            if i==4:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                ax.set_title(str(2004+i))
            elif i==7:
                ax.set_frame_on(False)
                ax.set_axis_off()
            else:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                    merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                ax.set_title(str(2004+i))

    return fig, axarr

def plot_stuff(xe, ye, H, cmap, grid, shelf_depth, ax):
    '''
    Do the main plotting stuff.
    '''

    XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))

    # Try with pcolor too
    mappable = ax.contourf(XE, YE, H.T, cmap=cmap, levels=np.linspace(0,100,11))
    ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

    return mappable


def plot_colorbar(fig, mappable, whichtype):
    '''
    Add colorbar to figure.
    '''

    # Make colorbar
    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')

    if whichtype == 'cross':
        cb.set_label('Probability of drifters crossing shelf (%)')
    elif 'coast' in whichtype: 
        cb.set_label('Probability of drifters reaching coastline region (%)')


def plot_finish(fig, whichtype, whichtime, shelf_depth):
    '''
    Save and close figure
    '''

    if whichtype == 'cross':
        fname = 'figures/' + whichtype + '/' + whichtime + str(shelf_depth) + '.png'
    elif 'coast' in whichtype: 
        fname = 'figures/' + whichtype + '/' + whichtime + '.png'

    fig.savefig(fname)
    plt.close(fig)

def plot_diff():
    '''
    Plot difference in transport between two histograms. Plot separately for convenience.
    '''

    shelf_depth = 20

    # Grid info
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, usebasemap=True)
    cmap = 'RdBu'

    d = np.load('figures/cross/seasonal20H.npz')
    Hboth = d['H']; xe=d['xe']; ye=d['ye']
    d.close()

    H = Hboth[0,:] - Hboth[1,:]

    fig = plt.figure(figsize=(6.8375, 6.6125))
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
    mappable = plot_stuff(xe, ye, H*100, cmap, grid, shelf_depth, ax)
    fig.colorbar(mappable)

    pdb.set_trace()



def run():

    # Which timing of plot: 'weatherband[1-3]', 'seasonal', 'interannual-winter', 'interannual-summer'
    whichtime = 'seasonal'#'interannual-winter'
    # Which type of plot: 'cross', 'coastCH', 'coastMX', 'coastLA', 
    #  'coastNTX', 'coastSTX', 'fsle', 'D2'
    whichtype = 'D2'

    shelf_depth = -20 # do 100 50 and 20 
    ishelf_depth = 0 # 2 1 0 index in cross array

    # Number of bins to use in histogram
    bins = (100,100) #(30,30)

    # Load in Files to read from based on which type of plot is being run
    Files, cmap = init(whichtime, whichtype)

    # Grid info
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, usebasemap=True)

    # Calculate xrange and yrange for histograms
    Xrange = [grid['xpsi'].min(), grid['xpsi'].max()]
    Yrange = [grid['ypsi'].min(), grid['ypsi'].max()]

    ## Calculate starting position histogram just once ##
    # Read in connectivity info (previously calculated). 
    # Drifters always start in the same place.
    pdb.set_trace()
    d = np.load(Files[0][0])
    # Histogram of starting locations
    if whichtype == 'cross': # results are in xg, yg
        xp, yp, _ = tracpy.tools.interpolate2d(d['xg0'], d['yg0'], grid, 'm_ij2xy')
    elif 'coast' in whichtype:  # results are in xp, yp
        xp = d['xp0']; yp = d['yp0']
    elif whichtype == 'D2': # results are in xg, yg
        # xp, yp are lonp, latp in this case
        xp, yp, _ = tracpy.tools.interpolate2d(d['xg'][:,0], d['yg'][:,0], grid, 'm_ij2ll')

    # For D2 and fsle, Hstart contains indices of drifters seeded in bins
    Hstart, xe, ye = calc_histogram(xp, yp, whichtype, bins=bins, Xrange=Xrange, Yrange=Yrange)

    d.close()

    # pdb.set_trace()

    # Set up overall plot
    fig, axarr = plot_setup(whichtime, grid) # depends on which plot we're doing

    # For D2 and fsle, H contains the metric calculation averaged over that bin
    H = np.zeros((len(Files), Hstart.shape[0], Hstart.shape[1]))

    # Loop through calculation files to calculate overall histograms
    # pdb.set_trace()
    for i, files in enumerate(Files): # Files has multiple entries, 1 for each subplot

        if whichtype == 'cross' or 'coast' in whichtype:
            Hcross = np.zeros(bins) # initialize
            #pdb.set_trace()
            HstartUse = Hstart*len(files) # multiply to account for each simulation


        for File in files: # now loop through the files for this subplot

            # Read in info
            d = np.load(File)
            if whichtype == 'cross': # results are in xg, yg
            # [number of depths,number of tracks] to store time of crossing or nan if it doesn't cross
                xg0 = d['xg0']; yg0 = d['yg0']
                cross = d['cross']
            elif 'coast' in whichtype:  # results are in xp, yp
                xp = d['xp0']; yp = d['yp0']
                conn = d['conn'] 
            elif whichtype == 'D2' or whichtype == 'fsle':
                xg = d['xg'][:]; yg = d['yg'][:]
            d.close()

            # Count the drifters for the shelf_depth that have a non-nan entry
            if whichtype == 'cross':
                ind = ~np.isnan(cross[ishelf_depth,:])
                xp, yp, _ = tracpy.tools.interpolate2d(xg0[ind], yg0[ind], grid, 'm_ij2xy')
            elif 'coast' in whichtype: 
                ind = ~np.isnan(conn)
                xp = xp[ind]; yp = yp[ind]
            elif whichtype == 'D2' or whichtype == 'fsle':
                lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')

            if whichtype == 'cross' or 'coast' in whichtype:
                # Calculate and accumulate histograms of starting locations of drifters that cross shelf
                Hcrosstemp, _, _ = calc_histogram(xp, yp, whichtype, bins=bins, Xrange=Xrange, Yrange=Yrange)
                Hcross = np.nansum( np.vstack((Hcross[np.newaxis,:,:], Hcrosstemp[np.newaxis,:,:])), axis=0)
            elif whichtype == 'D2' or whichtype == 'fsle':
                # Calculate the metric in each bin and combine for all files
                metric_temp, nnans = calc_metric(xp, yp, Hstart, whichtype)
                H[i,:] = H[i,:] + metric_temp
            # d.close()

        # Calculate overall histogram
        # pdb.set_trace()
        if whichtype == 'cross' or 'coast' in whichtype:
            H[i,:] = (Hcross/HstartUse)*100
        elif whichtype == 'D2':
            H[i,:] = H[i,:]/nnans
        elif whichtype == 'fsle':
            H[i,:] = 1./H[i,:]/nnans

        # Do subplot
        mappable = plot_stuff(xe, ye, H[i,:], cmap, grid, shelf_depth, axarr.flatten()[i])

    # save H
    if not os.path.exists('figures/' + whichtype): 
        os.makedirs('figures/' + whichtype)
    if not os.path.exists('figures/' + whichtype + '/' + whichtime): 
        os.makedirs('figures/' + whichtype + '/' + whichtime)

    if whichtype == 'cross':
        np.savez('figures/' + whichtype + '/' + whichtime + str(shelf_depth) + 'H.npz', H=H, xe=xe, ye=ye)
    elif 'coast' in whichtype: 
        np.savez('figures/' + whichtype + '/' + whichtime + 'H.npz', H=H, xe=xe, ye=ye)

    # Add colorbar
    plot_colorbar(fig, mappable, whichtype)
    # pdb.set_trace()

    # save and close
    plot_finish(fig, whichtype, whichtime, shelf_depth)


if __name__ == "__main__":
    run()    
