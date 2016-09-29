"""
Use indices from shelfconn cross shelf transport analysis files to count
how many of the drifters also exit the domain
"""

import numpy as np
import netCDF4 as netCDF
from glob import glob
import os

Files = glob('calcs/shelfconn/*gc.npz')

shelf_depths = [20, 50, 100, 200, 300, 400, 500]
idepth = 2

nnansmean = 0
for i, calcfilename in enumerate(Files):
    # print File

    # file where tracks are stored
    name = calcfilename.split('/')[2]

    # file where we will store the drifter exit info
    crossexitfilename = 'calcs/shelfconn/exitdomain/' + name

    print crossexitfilename
    # don't repeat calc if last file was created
    if os.path.exists(crossexitfilename):
        # continue
        print 'adding to mean'
        d = np.load(crossexitfilename)
        nnansmean += d['nnans']
        if i==0:
            t = d['t']
        d.close()

    else:  # do analysis for the simulation
        # file where cross-shelf info is stored
        calcfile = np.load(calcfilename)

        trackfilename = 'tracks/' + name[:-3] + 'nc'
        trackfile = netCDF.Dataset(trackfilename)

        # indices of drifters that cross the shelf from this file
        inds = ~np.isnan(calcfile['cross'][idepth, :])

        # load in drifter tracks
        # assume we only need xg and not yg also
        xg = trackfile.variables['xg'][:]
        xg = xg[inds]
        irem = xg == -1
        xg[irem] = np.nan

        # want number of active drifters in time
        nnans = (~np.isnan(xg)).sum(axis=0)
        t = trackfile.variables['tp'][:]

        # save exit drifter info
        np.savez(crossexitfilename, nnans=nnans, t=t)
        trackfile.close()
        calcfile.close()

nnansmean /= len(Files)
np.savez('calcs/shelfconn/exitdomain/overallaverage.npz', nnansmean=nnansmean, t=t)