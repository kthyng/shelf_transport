'''
Script to run drifters at 1km initial spacing daily forward for 30 days.
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
from tracpy.tracpy_class import Tracpy


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
# currents_filename = list(np.sort(glob.glob('/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_????.nc')))

# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# can't aggregate years between 2012 and before with 2013 and 2014 bc they have different variables
# years = np.arange(2011,2013)
# currents_filename = []
# for year in years:
#     currents_filename.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_????.nc')))

years = np.arange(2004,2005)
currents_filename = []
for year in years:
    currents_filename.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_*.nc')))

grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)


def init(name):
    '''
    Initialization for the simulation.
    '''

    # loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'

    time_units = 'seconds since 1970-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 25 # in-between tracks: 12 # old tracks: 25 

    # Number of steps to divide model output for outputting drifter location
    N = 5

    # Number of days
    ndays = 30

    # This is a forward-moving simulation
    ff = 1 

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0. # old tracks: 5.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines.
    dostream = 0

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid_filename=grid_filename, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, dostream=dostream, savell=False, doperiodic=0, 
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, 
                time_units=time_units, usebasemap=True, grid=grid, vert_filename=vert_filename)

    # tp._readgrid()

    # initial separation distance of drifters, in meters, from sensitivity project
    dx = 500
    seedsfile = 'calcs/seeds' + str(dx) + '.npz'
    if os.path.exists(seedsfile):
        seeds = np.load(seedsfile)
        lon0 = seeds['lon0']; lat0 = seeds['lat0']
        seeds.close()
    else:
        # Initial lon/lat locations for drifters
        # Start uniform array of drifters across domain using x,y coords
        llcrnrlon = tp.grid['lonr'].min(); urcrnrlon = tp.grid['lonr'].max(); 
        llcrnrlat = tp.grid['latr'].min(); urcrnrlat = tp.grid['latr'].max(); 
        xcrnrs, ycrnrs = tp.grid['basemap']([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
        X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], dx), np.arange(ycrnrs[0], ycrnrs[1], dx))

        lon0, lat0 = tp.grid['basemap'](X, Y, inverse=True)

        # Eliminate points that are outside domain or in masked areas
        lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)

        # save starting locations for future use
        np.savez(seedsfile, lon0=lon0, lat0=lat0)

    # # equal weightings for drifters for transport.
    # T0 = np.ones(lon0.size, order='F')

    # U = np.ma.zeros(tp.grid['xu'].shape, order='F')
    # V = np.ma.zeros(tp.grid['xv'].shape, order='F')

    # pdb.set_trace()

    return tp, lon0, lat0


def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    overallstartdate = datetime(2004, 1, 1, 0, 1)
    overallstopdate = datetime(2004, 1, 2, 0, 1)
    # overallstopdate = datetime(2014, 7, 1, 4, 1)

    date = overallstartdate

    # Start from the beginning and add days on for loop
    # keep running until we hit the next month
    while date < overallstopdate:

        name = date.isoformat()[0:13]

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
            not os.path.exists('tracks/' + name + 'gc.nc'):

            # Read in simulation initialization
            tp, lon0, lat0 = init(name)

            # Run tracpy
            # Save directly to grid coordinates
            lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

        # Increment by 24 hours for next loop, to move through more quickly
        date = date + timedelta(hours=24)


if __name__ == "__main__":
    run()    
