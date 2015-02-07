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
import op
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar



# of the drifters that at some point enter the coastal boxes
# When do they first enter a box? Integrate together for likelihood map
# How many drifters enter each box and at what original time? Save original location


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

# load in initial drifter starting locations in grid space
d = np.load('calcs/xyp0.npz')
xp0 = d['xp0']; yp0 = d['yp0']
d.close()

bins = (100,100)
# Calculate xrange and yrange for histograms
Xrange = [grid['xpsi'].min(), grid['xpsi'].max()]
Yrange = [grid['ypsi'].min(), grid['ypsi'].max()]

# Save a histogram of the number of drifters that started in each bin
Hstart, xe, ye = np.histogram2d(xp0, yp0, bins=bins, 
                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])

# # Find indices of all drifters that start in the coastal boxes
# # start indices of drifters that start in each box path
pts = np.load('calcs/alongcoastconn/inds-in-coast-paths.npz')['pts']
# pt = [] # aggregated indices of drifters that start in coast boxes
# [pt.extend(pts[j]) for j in xrange(len(pts))]
# xp0coast = xp0[pt]; yp0coast = yp0[pt]


def likelihood():
    '''
    Aggregate likelihood of connection from locations with coast 
    in different time periods.
    '''

    nd = np.load('calcs/xyg0.npz')['xg0'].size # # of drifters

    # Loop through along-coast boxes to find which other boxes they are connected to
    years = np.arange(2004,2015)
    months = [1,2,7,8]
    days = np.array([3,5,10,15,20,30])
    for year in years:
        for month in months:
            fname = 'calcs/coastconn/likelihood/hist-' + str(year) + '-' + str(month).zfill(2) + '.npz'
            if not os.path.exists(fname):
                Files = glob.glob('calcs/alongcoastconn/' + str(year) \
                            + '-' + str(month).zfill(2) + '-*T0*.npz')

                # # likelihood histogram for each advection time examined
                # H = np.zeros((days.size, bins[0], bins[1]))
                # number of drifters reaching each coast box (yes or no, no weighting) in day days
                ndbox = np.zeros((days.size, len(pts)))
                # Histograms of where drifters originate for advection times by coast box
                H = np.zeros((days.size, len(pts), bins[0], bins[1]))
                
                for File in Files:
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

                    # # necesary to translate between points in different arrays and have array size correct
                    # xp0temp = xp0[code][np.newaxis,:]
                    # xp0temp = xp0temp.repeat(inbox.shape[0], axis=0)
                    # yp0temp = yp0[code][np.newaxis,:]
                    # yp0temp = yp0temp.repeat(inbox.shape[0], axis=0)

                    for i, day in enumerate(days): # loop through number of advection days

                        # Drifters that enter a coast box within day days [coast box x set of drifters]
                        ind = (inbox[:,:,0]<=day)

                        # How many drifters enter each box by day days?
                        ndbox[i,:] += ind.sum(axis=1)

                        # loop through each coast box path to calculate seperate origin histogram
                        for j in xrange(len(pts)):

                            # indices of drifters that start in this box, referenced to shelf transport seeding
                            pt = pts[j] 

                            # projected drifter origin locations that end up in this coast box in day days
                            xpsave = xp0[code][ind[j,:]]; ypsave = yp0[code][ind[j,:]]

                            # # Original xp, yp locations of drifters that enter a coast box within day days
                            # xpsave = xp0temp[ind]; ypsave = yp0temp[ind]
                            # Add in drifters that start in coast boxes
                            xpsave = np.concatenate((xpsave, xp0[pt]))
                            ypsave = np.concatenate((ypsave, yp0[pt]))

                            # Find histogram of xpsave, ypsave points for this simulation/box origin points
                            Htemp, _, _ = np.histogram2d(xpsave, ypsave, bins=bins, 
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                            # aggregate number of drifters starting in different histogram bins 
                            # that reach coastline for each month/year combination
                            H[i,j,:,:] += Htemp 



                # Save the month/year's worth of histograms
                # numfiles is to calculate the number drifters from bins for the the number of runs
                # aggregated together, compared with the appropriate number of starting drifters overall
                np.savez(fname, H=H, xe=xe, ye=ye, days=days, numfiles=len(Files), ndbox=ndbox)
                # pdb.set_trace()
            
        

# if __name__ == "__main__":
#     run()     
