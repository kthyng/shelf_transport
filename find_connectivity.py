'''
Script to do analysis for what drifters cross the shelf. This will be run on hafen,
where all of the tracks files are.

Run from shelf_transport directory.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob
from matplotlib.mlab import find

# Find files to run through
Files = glob.glob('/Volumes/Emmons/projects/gisr/tracks/all_f/*gc.npz')
# Files = glob.glob('tracks/*gc.nc')

shelf_depth = 100

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)

# Loop through files for analysis
for File in Files:

    if os.path.exists('calcs/shelfconn/' + File.split('/')[-1][:-5] + '.npz'):
        continue

    # Get necessary info from File
    d = np.load(File)
    xg = d['xg']; yg = d['yg']; tg = d['tg']
    # d = netCDF.Dataset(File)
    # xg = d.variables['xg'][:]
    # yg = d.variables['yg'][:]
    # tg = d.variables['tg'][:]
    d.close()

    # Time information
    days = (tg-tg[0])/(3600.*24)
    iday = find(days==2) - find(days==1) # number of indices per day

    # Calculate the depths for the drifter positions
    nanind = np.isnan(xg) # indices where nans are location in xg, yg; for reinstitution of nans
    xg[nanind] = 0; yg[nanind] = 0 # need to preserve the shape of xg, yg, and have valid indices to use with grid['h']
    h = grid['h'][np.ceil(xg).astype(int), np.ceil(yg).astype(int)] # depths at all drifter locations
    h[nanind] = np.nan # reinstitute nans into depths of trajectories

    # Change to projected drifter locations now
    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
    xp[nanind] = np.nan; yp[nanind] = np.nan
    del(xg,yg) # don't need grid info anymore

    # Initial array of information
    ntrac = xp.shape[0] # number of drifters
    ndays = np.arange(1,31) # number of simulation days
    cross = np.ones((ntrac, ndays.size))*np.nan # to store analysis

    # Loop through number of days
    for nday in ndays:

        nind = find(nday>=days)[-1] # number of indices into time dimension this nday corresponds to

        # find indices of drifters deeper and shallower than shelf_depth)
        ideep = h[:,:nind] > shelf_depth # indices of deeper locations
        ishallow = h[:,:nind] <= shelf_depth # indices of shallower locations

        # Find consecutive values of array (both shallower and deeper than shelf_depth). Searches in time direction
        # from http://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy
        pdb.set_trace()
        adeep = np.array_split(ideep, np.where(np.diff(ideep)!=1, axis=1)[0]+1, axis=1) # gives back list of sets of indices of consecutive deeper locations
        ashallow = np.array_split(ishallow, np.where(np.diff(ishallow)!=1, axis=1)[0]+1, axis=1) # gives back list of sets of indices of consecutive deeper locations

        # Search deep and shallow info to see if on both sides for enough time


        # Using drifter locations in grid coordinates, make drifter-sized array of 1's and 0's, 1 if the depth is both above and below
        # the shelf_depth during the nday length of considered simulation, and 0 if not. This is a function of nday and shelf_depth.
        # crosses = np.zeros(T0.shape) # elements will be changed to 1 if a given drifter is both above and below isobath at some point
        hsum = np.sum(h[:,:nind]>shelf_depth, axis=1) # ONLY LOOK UP THROUGH NDAY IN TIME FOR CROSSINGS
        # need to both have at least one point greater than shelf_depth but not all non-nan positions have over shelf_depth depth
        crosses = (hsum>0)*(hsum<np.sum(~hind[:,:nind], axis=1)) # bool of whether drifter crosses shelf



