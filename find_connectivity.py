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
Files = glob.glob('/Volumes/Emmons/projects/gisr/tracks/all_f/*gc.npz')[0:2]
# Files = glob.glob('tracks/*gc.nc')

shelf_depth = 100

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)

# Loop through files for analysis
for File in Files:

    savefile = 'calcs/shelfconn/' + File.split('/')[-1][:-5] + '.npz'

    if os.path.exists(savefile):
        continue

    # Get necessary info from File
    d = np.load(File)
    xg = d['xg']; yg = d['yg']; tg = d['tg']
    # d = netCDF.Dataset(File)
    # xg = d.variables['xg'][:]
    # yg = d.variables['yg'][:]
    # tg = d.variables['tp'][:]
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
    ndays = np.arange(0,31) # number of simulation days # the first few days will be empty
    cross = np.zeros((ntrac, ndays.size)) # to store analysis. True if crosses shelf by nday.

    # Loop through drifters
    for i in xrange(ntrac):

        # find indices of drifters deeper and shallower than shelf_depth)
        ideep = h > shelf_depth # indices of deeper locations
        ishallow = h <= shelf_depth # indices of shallower locations
        # pdb.set_trace()

        # Loop through number of days
        for nday in ndays:

            # nothing will happen these days
            if nday==0:
                continue
            if nday==1:
                continue

            nind1 = find(nday>=days)[-1] # number of indices into time dimension this nday corresponds to
            nind2 = find(~np.isnan(h[i,:]))[-1] # last non-nan
            # want the last index before nan's start or the appropriate timing for nday, whichever comes first
            nind = min((nind1,nind2))

            # Find consecutive values of array (both shallower and deeper than shelf_depth). Searches in time direction
            # from http://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy
            # number of entries in adeep list is the number of times it changes sides of the shelf, and number of entries
            # in each entry of list is the number of outputs it is there. True is deep and False is shallow.
            # nan's show up as False
            adeep = np.array_split(ideep[i,:nind], np.where(np.diff(ideep[i,:nind])==1)[0]+1)
            ashallow = np.array_split(ishallow[i,:nind], np.where(np.diff(ishallow[i,:nind])==1)[0]+1)

            # Search through chunks of consecutive true's and false's (deep/shallow) to see if on both sides for enough time
            # move on if find a true
            for j in xrange( len(adeep) ):
            # for j in xrange( max(len(adeep),len(ashallow)) ):

                # Find when starts deep and has more than iday elements on shallow side
                if (h[i,0]>shelf_depth) and (ashallow[j].astype(int).sum()>iday):
                    cross[i,nday] = True
                    continue

                # Find when starts shallow and has more than iday elements on deep side
                elif (h[i,0]<shelf_depth) and (adeep[j].astype(int).sum()>iday):
                    cross[i,nday] = True
                    continue

            # If it has crossed for nday, then it will cross for the rest of them too
            if cross[i,nday] == True:
                cross[i,nday+1:] = True
                continue
            # If the track has nan'ed out at this point, whatever the cross value is will be true for the rest of ndays
            elif nind2 < nind1:
                cross[i,nday+1:] = cross[i,nday]
                continue

    # save array
    np.savez(savefile, cross=cross)


