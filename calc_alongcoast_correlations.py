
'''
Calculate correlations between alongcoast connectivity (seasonal, interannual)
and forcing mechanisms along the coast.
'''

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
# import tracpy
from matplotlib.path import Path
import cartopy
import os
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import pandas as pd
from glob import glob


## read in alongcoast connectivity line (summed across matrix) ##
# want seasonal and interannual
# work from jupyter notebook in shelf_transport/notebooks/porta_tofrom_surfside.ipynb
base = 'calcs/alongcoastconn/'
lname = base + 'lines.npz'
d = np.load(lname)
# gives upcoast alongcoast conn as function of starting position
startingp = d['startingp'];
# gives downcoast alongcoast conn as function of starting position
startingn = d['startingn']
endingp = d['endingp']; endingn = d['endingn']; dates = d['dates']; dist = d['dist']
d.close()
# startingp, etc, are accidentally too long
sp = pd.DataFrame(index=[pd.Timestamp(date) for date in dates], columns=np.arange(0,len(dist)), data=startingp[:dates.size])
sn = pd.DataFrame(index=[pd.Timestamp(date) for date in dates], columns=np.arange(0,len(dist)), data=startingn[:dates.size])
ep = pd.DataFrame(index=[pd.Timestamp(date) for date in dates], columns=np.arange(0,len(dist)), data=endingp[:dates.size])
en = pd.DataFrame(index=[pd.Timestamp(date) for date in dates], columns=np.arange(0,len(dist)), data=endingn[:dates.size])


# read in connection between PA and SS
# work from jupyter notebook in shelf_transport/notebooks/porta_tofrom_surfside.ipynb
base = 'calcs/alongcoastconn/conn_in_time/pa_ss/'
t = np.load(base + 't.npz')['t']  # in decimal days
oname = base + 'aggregated.npz'
Files = glob(base + '2*.npz')
dates = [pd.Timestamp(File.split('/')[-1][:-4]) for File in Files]
d = np.load(oname)
# to save results, [start datetime (dates), 30 days every 4 hours]
ss2pa = d['ss2pa']
pa2ss = d['pa2ss']
d.close()


## read in alongcoast geometric vectors ##
# from shelf_transport/calc_alongcoast_windangle.py
vecname = 'calcs/along_coast_wind/coast_vectors.npz'
d = np.load(vecname)
# upcoast-pointing vector, goes from MX to LA. + is upcoast and toward land.
veccoast = d['veccoast']
# vector that is rotated 90 degrees off from coast direction
veccoast90 = d['veccoast90']
# along-coast distance (km)
dist = d['dist']
# xgcoast=xgcoast, ygcoast=ygcoast, loncoast=loncoast, latcoast=latcoast)


## Select model grid points that are within the coastal boxes ##
inds_rho_boxes = np.load('calcs/coastpaths_pts.npz')['inds_rho_boxes']
inds_u_boxes = np.load('calcs/coastpaths_pts.npz')['inds_u_boxes']
inds_v_boxes = np.load('calcs/coastpaths_pts.npz')['inds_v_boxes']


## read in along and across coast wind component and river water (seasonal and interannual)
# calculated in calc_alongcoast_wind_river.py
base_wind = 'calcs/along_coast_wind/'
base_river = 'calcs/along_coast_rivers/'

years = np.arange(2004, 2015)
year = 2004
# for year in years:
    # dates = m['ocean_time'].sel(ocean_time=str(year))
fname_wind = base_wind + str(year) + '.csv'
fname_river = base_river + str(year) + '.csv'

walong = pd.read_csv(fname_wind, parse_dates=True, index_col=0)
wacross = pd.read_csv(fname_wind, parse_dates=True, index_col=0)
miss = pd.read_csv(fname_river, parse_dates=True, index_col=0)
atch = pd.read_csv(fname_wind, parse_dates=True, index_col=0)
braz = pd.read_csv(fname_wind, parse_dates=True, index_col=0)



## Look at spatial correlations
