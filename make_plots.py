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

    base = 'calcs/shelfconn/'
    #pdb.set_trace()
    if whichtime=='seasonal':
        Files = []
        # Seasonal returns Files that has two entries of lists
        Files.append(glob(base + '*-0[6-7]-*.npz'))
        #Files.extend(glob(base + '*-08-*.npz'))
        Files.append(glob(base + '*-0[1-2]-*.npz'))
        #Files.extend(glob(base + '*-02-*.npz'))

    if whichtype=='cross':
        cmap = 'YlOrRd'

    return Files, cmap


def calc_histogram(xp, yp, bins=(60,60), 
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

    if Xrange!=None and Yrange!=None:
        H, xe, ye = np.histogram2d(xp, yp, bins=bins, 
                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
    else:
        H, xe, ye = np.histogram2d(xp, yp, bins=bins)

    return H, xe, ye


def plot_setup(whichtime, grid):
    '''
    Set up the figure for plotting based on the string in which.

    Inputs:

    whichtime           Which type of plot to do
    grid            Previously-read-in grid dictionary, with basemap included
    '''

    if whichtime=='seasonal':
        # fig = plt.figure(figsize=(22,10))
        fig, axarr = plt.subplots(2)#, sharex=True)
        # ax = fig.add_subplot(1,2,1)
        for ax in axarr:
            tracpy.plotting.background(grid=grid, ax=ax)
            # Titles for subplots
            ax.set_title('Winter')
            ax.set_title('Summer')
            # suptitle
            fig.suptitle('Probability of material crossing the shelf in 30 days, 2004-2010', y=.94)


    return fig, axarr

def plot_stuff(xe, ye, H, cmap, grid, shelf_depth, ax):
    '''
    Do the main plotting stuff.
    '''

    XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))

    # Try with pcolor too
    ax.contourf(XE, YE, H.T, cmap, levels=np.linspace(0,100,11))
    ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)


def plot_colorbar(fig):
    '''
    Add colorbar to figure.
    '''

    # Make colorbar
    # Horizontal colorbar below plot
    cax = fig.add_axes([0.25, 0.05, 0.48, 0.02]) #colorbar axes
    cb = colorbar(cax=cax,orientation='horizontal')
    cb.set_label('Probability (%)')
    d.close()


def plot_finish(fig, whichtype, whichtime):
    '''
    Save and close figure
    '''

    fname = 'figures/' + whichtype + '/' + whichtime + '.png'

    if not os.path.exists('figures/' + whichtype): 
        os.makedirs('figures/' + whichtype)
    if not os.path.exists('figures/' + whichtype + '/' + whichtime): 
        os.makedirs('figures/' + whichtype + '/' + whichtime)

    fig.savefig(fname)
    plt.close(fig)

def run():

    # Which timing of plot: 'weatherband', 'seasonal', 'interannual'
    whichtime = 'seasonal'
    # Which type of plot: 'cross' or 'coast'
    whichtype = 'cross'

    shelf_depth = 100 # do 50 also
    ishelf_depth = 2 # index in cross array

    # Number of bins to use in histogram
    bins = (60,60)

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
    d = np.load(Files[0][0])
    # Histogram of starting locations
    Hstart, xe, ye = calc_histogram(d['xg0'], d['yg0'], bins=bins, Xrange=Xrange, Yrange=Yrange)
    d.close()

    # Set up overall plot
    fig, axarr = plot_setup(whichtime, grid) # depends on which plot we're doing

    # Loop through calculation files to calculate overall histograms
    # pdb.set_trace()
    for i, files in enumerate(Files): # Files has multiple entries, 1 for each subplot

        Hcross = np.zeros(bins) # initialize
        #pdb.set_trace()
        Hstart *= len(files) # multiply to account for each simulation

        for File in files: # now loop through the files for this subplot

            # Read in connectivity info (previously calculated)
            d = np.load(File)
            xg0 = d['xg0']; yg0 = d['yg0']
            # [number of depths,number of tracks] to store time of crossing or nan if it doesn't cross
            cross = d['cross']
            d.close()

            # Count the drifters for the shelf_depth that have a non-nan entri
            ind = ~np.isnan(cross[ishelf_depth,:])

            # Calculate and accumulate histograms of starting locations of drifters that cross shelf
            Hcrosstemp, _, _ = calc_histogram(xg0[ind], yg0[ind], bins=bins, Xrange=Xrange, Yrange=Yrange)
            Hcross = np.nansum( np.vstack((Hcross, Hcrosstemp)))

        # Calculate overall histogram
        H = Hstart/Hcross

        # Do subplot
        pdb.set_trace()
        plot_stuff(xe, ye, H, cmap, grid, shelf_depth, axarr[i])

    # Add colorbar
    plot_colorbar(fig)

    # save and close
    plot_finish(fig, whichtime, whichtype)


if __name__ == "__main__":
    run()    
