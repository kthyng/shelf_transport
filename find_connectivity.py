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

# 'cross' or 'coast'
which = 'cross'

def load_tracks(File):
    '''
    Load in and return tracks from File
    '''

    # Get necessary info from File
    d = netCDF.Dataset(File)
    xg = d.variables['xg'][:]
    yg = d.variables['yg'][:]
    tp = d.variables['tp'][:]
    d.close()

    return xg, yg, tp



def run():


    # Find files to run through
    # Files = glob.glob('/Volumes/Emmons/projects/gisr/tracks/all_f/*gc.npz')[0:2]
    Files = glob.glob('tracks/*gc.nc')

    shelf_depths = [20, 50, 100, 200, 300, 400, 500]

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc)

    # Loop through files for analysis
    for File in Files:

        if which=='cross':

            shelffile = 'calcs/shelfconn/' + File.split('/')[-1][:-3] + '.npz'

            if os.path.exists(shelffile): # don't repeat calc if last file was created
                continue

            xg, yg, tp = load_tracks(File)

            # Time information
            days = (tp-tp[0])/(3600.*24)
            iday = find(days==2) - find(days==1) # number of indices per day
            iday = int(iday) 

            # Calculate the depths for the drifter positions
            nanind = np.isnan(xg)# + xg==-1 # indices where nans are location in xg, yg; for reinstitution of nans
            xg[nanind] = 0; yg[nanind] = 0 # need to preserve the shape of xg, yg, and have valid indices to use with grid['h']
            #pdb.set_trace()
            h = grid['h'][np.ceil(xg).astype(int), np.ceil(yg).astype(int)] # depths at all drifter locations
            h[nanind] = np.nan # reinstitute nans into depths of trajectories

            # Initial array of information
            ntrac = xg.shape[0] # number of drifters
            ndays = np.arange(0,31) # number of simulation days # the first few days will be empty
            cross = np.ones((len(shelf_depths), ntrac))*np.nan # [number of depths,number of tracks] to store time of crossing or nan if it doesn't cross

            # Cross will hold timing of when the drifter crosses the shelf the first time if it 
            # spends at least iday points both deeper and shallower than shelf_depth, and nan otherwise
            for i,shelf_depth in enumerate(shelf_depths):
                ideep = h>shelf_depth
                ndeep = np.sum(ideep, axis=1) # number deeper than shelf depth
                ishallow = h<shelf_depth
                nshallow = np.sum(ishallow, axis=1) # number shallower than shelf depth
                # need to have at least iday points deeper and shallower than shelf_depth 
                icross = (ndeep>iday) * (nshallow>iday) # bools, True if crosses
                # pdb.set_trace()
                whencross = np.diff(ideep, axis=1) # will pick out when the drifters change from deep to shallow or vice versa
                # days2d = np.expand_dims(days,axis=0).repeat(ntrac, axis=0)
                # picks out the first time the switch happens, and gives the time in days
                # is zero if not
                iwhen = (days[1:]*whencross).argmax(axis=1) 
                # when = days[(days[1:]*whencross).argmax(axis=1)+1]
                # contains the time in days when the shelf is crossed the first time
                # has to be on both sides of the shelf for idays
                # is nan if doesn't cross the shelf
                cross[i,icross] = days[iwhen[icross]] 

            # save array
            np.savez(shelffile, xg0=xg[:,0], yg0=yg[:,0], cross=cross)


        elif which=='coast':

            MXfile = 'calcs/coastconn/MX/' + File.split('/')[-1][:-4] + '.npz'
            STXfile = 'calcs/coastconn/STX/' + File.split('/')[-1][:-4] + '.npz'
            NTXfile = 'calcs/coastconn/NTX/' + File.split('/')[-1][:-4] + '.npz'
            CHfile = 'calcs/coastconn/CH/' + File.split('/')[-1][:-4] + '.npz'
            LAfile = 'calcs/coastconn/LA/' + File.split('/')[-1][:-4] + '.npz'
            #pdb.set_trace()
            if os.path.exists(LAfile): # don't repeat calc if last file was created
                continue

            xg, yg, tp = load_tracks(File)
 
            # Time information
            days = (tp-tp[0])/(3600.*24)
            iday = find(days==2) - find(days==1) # number of indices per day
            iday = int(iday) 
            # pdb.set_trace()
            # Load in coastline region info
            # Mexico
            dconn = np.load('calcs/MXpts.npz')
            xp = np.asarray(dconn['MX'])[:,0]; yp = np.asarray(dconn['MX'])[:,1]
            MXpath = Path(np.vstack((xp, yp)).T)
            dconn.close()
            # S Texas
            dconn = np.load('calcs/STXpts.npz')
            xp = np.asarray(dconn['STX'])[:,0]; yp = np.asarray(dconn['STX'])[:,1]
            STXpath = Path(np.vstack((xp, yp)).T)
            dconn.close()
            # N Texas
            dconn = np.load('calcs/NTXpts.npz')
            xp = np.asarray(dconn['NTX'])[:,0]; yp = np.asarray(dconn['NTX'])[:,1]
            NTXpath = Path(np.vstack((xp, yp)).T)
            dconn.close()
            # Chenier
            dconn = np.load('calcs/CHpts.npz')
            xp = np.asarray(dconn['CH'])[:,0]; yp = np.asarray(dconn['CH'])[:,1]
            CHpath = Path(np.vstack((xp, yp)).T)
            dconn.close()
            # Louisiana
            dconn = np.load('calcs/LApts.npz')
            xp = np.asarray(dconn['LA'])[:,0]; yp = np.asarray(dconn['LA'])[:,1]
            LApath = Path(np.vstack((xp, yp)).T)
            dconn.close()
            # dconn = np.load('calcs/LAinside.npz')
            # LA = dconn['inside'] # indices for flattened rho array
            # dconn.close()

            # Change to projected drifter locations now
            nanind = np.isnan(xg) #+ xg==-1 # indices where nans are location in xg, yg; for reinstitution of nans
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
            xp[nanind] = np.nan; yp[nanind] = np.nan
            del(xg,yg) # don't need grid info anymore

            # Initial array of information
            # cross = np.zeros((ntrac, ndays.size), dtype=bool) # to store analysis. True if crosses shelf by nday.
            ntrac = xp.shape[0]            
            MX = np.ones(ntrac)*np.nan # to store analysis. 
            STX = np.ones(ntrac)*np.nan # to store analysis. 
            NTX = np.ones(ntrac)*np.nan # to store analysis. 
            CH = np.ones(ntrac)*np.nan # to store analysis. 
            LA = np.ones(ntrac)*np.nan # to store analysis. 

            # which points are inside the region
            # Mexico
            #pdb.set_trace()
            inside = MXpath.contains_points(np.vstack((xp.flat, yp.flat)).T).reshape(xp.shape)
            iinside = inside.sum(axis=1).astype(bool) # true if goes inside
            whencross = np.diff(inside, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            iwhen = (days[1:]*whencross).argmax(axis=1) 
            MX[iinside] = days[iwhen[iinside]]

            # S Texas
            inside = STXpath.contains_points(np.vstack((xp.flat, yp.flat)).T).reshape(xp.shape)
            iinside = inside.sum(axis=1).astype(bool) # true if goes inside
            whencross = np.diff(inside, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            iwhen = (days[1:]*whencross).argmax(axis=1) 
            STX[iinside] = days[iwhen[iinside]]

            # N Texas
            inside = NTXpath.contains_points(np.vstack((xp.flat, yp.flat)).T).reshape(xp.shape)
            iinside = inside.sum(axis=1).astype(bool) # true if goes inside
            whencross = np.diff(inside, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            iwhen = (days[1:]*whencross).argmax(axis=1) 
            NTX[iinside] = days[iwhen[iinside]]

            # Chenier
            inside = CHpath.contains_points(np.vstack((xp.flat, yp.flat)).T).reshape(xp.shape)
            iinside = inside.sum(axis=1).astype(bool) # true if goes inside
            whencross = np.diff(inside, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            iwhen = (days[1:]*whencross).argmax(axis=1) 
            CH[iinside] = days[iwhen[iinside]]

            # Louisiana
            inside = LApath.contains_points(np.vstack((xp.flat, yp.flat)).T).reshape(xp.shape)
            iinside = inside.sum(axis=1).astype(bool) # true if goes inside
            whencross = np.diff(inside, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            iwhen = (days[1:]*whencross).argmax(axis=1) 
            LA[iinside] = days[iwhen[iinside]]

            # save array
            np.savez(MXfile, xp0=xp[:,0], yp0=yp[:,0], conn=MX)
            np.savez(STXfile, xp0=xp[:,0], yp0=yp[:,0], conn=STX)
            np.savez(NTXfile, xp0=xp[:,0], yp0=yp[:,0], conn=NTX)
            np.savez(CHfile, xp0=xp[:,0], yp0=yp[:,0], conn=CH)
            np.savez(LAfile, xp0=xp[:,0], yp0=yp[:,0], conn=LA)


if __name__ == "__main__":
    run()    
