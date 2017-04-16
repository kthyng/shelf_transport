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
from matplotlib.mlab import find, Path
import time

def load_tracks(File):
    '''
    Load in and return tracks from File
    '''

    # Get necessary info from File
    d = netCDF.Dataset(File)
    xg = d.variables['xg'][:]
    yg = d.variables['yg'][:]
    #tp = d.variables['tp'][:]
    # Calculate times since they are messed up if a drifter exits the domain
    tseas = d.variables['tseas'][:]
    N = d.variables['N'][:]
    iday = int((3600.*24)/(tseas/N))
    days = np.linspace(0, 30, xg.shape[1]) #tseas/N/(3600.*24)) # forward or backward, but doesn't matter here I guess
    # pdb.set_trace()
    d.close()

    return xg, yg, iday, days



def run():

    # Find files to run through
    # Files = glob.glob('tracks/2014-0[1,2]-*gc.nc')
    # Files.extend(glob.glob('tracks/2014-0[7,8]-*gc.nc'))
    Files = glob.glob('tracks/2008-??-*gc.nc')

    # grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    # vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
    # grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

    # load in paths for coastline boxes
    d = np.load('calcs/coastpaths.npz') # use paths in grid space
    paths = d['pathsg']
    pathouter = d['outerpathg'].item()
    d.close()

    # Loop through files for analysis
    for File in Files:

        savefile = 'calcs/alongcoastconn/' + File.split('/')[-1][:-5] + '.npz'
        if os.path.exists(savefile): # don't repeat calc if last file was created
            continue

        # import pdb; pdb.set_trace()
        # print(File)
        # import time
        # stime = time.time()
        xg, yg, iday, days = load_tracks(File)

        # Change to projected drifter locations now
        nanind = (xg==-1) # indices where nans are location in xg, yg; for reinstitution of nans
        xg[nanind] = np.nan; yg[nanind] = np.nan

        # Remove points that never enter the outer box paths
        inouter = pathouter.contains_points(np.vstack((xg.flat, yg.flat)).T).reshape(xg.shape)
        # print('made it past contains_points')
        # print(time.time()-stime)
        iinside = find(inouter.sum(axis=1).astype(bool)) # indices of drifters that go inside
        # print '0 time for outer path comparison: ', time.time() - tic_start
        # save a pointer of these drifters within the large arrays of drifters
        xgin = xg[iinside,:]; ygin = yg[iinside,:]

        # Initial array of information
        ncrosses = 5
        inbox = np.ones((len(paths), iinside.size, ncrosses))*np.nan # to store analysis.
        outbox = np.ones((len(paths), iinside.size, ncrosses))*np.nan # to store analysis.

        for i,path in enumerate(paths):

            # which points are inside the regions
            inside = path.contains_points(np.vstack((xgin.flat, ygin.flat)).T).reshape(xgin.shape)
            whencross = np.diff(inside.astype(int), axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            whencross_sorted = np.sort(whencross, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            isort = np.argsort(whencross, axis=1)
            outbox[i,:,:] = days[isort[:,:ncrosses]] # allow for up to ncrosses re-exits
            # nan out the entries that aren't actually exiting
            iout = whencross_sorted[:,:ncrosses]!=-1
            outbox[i,iout] = np.nan
            outbox[i,:,:] = np.sort(outbox[i,:,:], axis=1)

            inbox[i,:,:] = days[isort[:,-ncrosses:]] # allow for up to ncrosses re-enters
            # nan out the exits that aren't actually exiting
            iin = whencross_sorted[:,-ncrosses:]!=1
            inbox[i,iin] = np.nan
            inbox[i,:,:] = np.sort(inbox[i,:,:], axis=1)

        # save array
        np.savez(savefile, inbox=inbox, outbox=outbox, iinside=iinside)



if __name__ == "__main__":
    run()
