'''
File to make table for multiple variable correlation with shelf transport.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
from datetime import datetime, timedelta
import glob
from cmPong import cmPong
from matplotlib.mlab import find
import bisect
from matplotlib import delaunay
import csv


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


# Loop through years
years = np.arange(2004,2011)


for year in years:

    ## River forcing ##
    r = netCDF.Dataset('/atch/raid1/zhangxq/Projects/txla_nesting6/TXLA_river_4dyes_2011.nc')
    # River timing
    unitsRiver = r.variables['river_time'].units
    datesRiver = netCDF.num2date(r.variables['river_time'][:], unitsRiver)
    tRiver = r.variables['river_time'][:]
    # all of river input
    Q = np.abs(r.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
    # # start and end dates for river discharge
    # tstartRiver = netCDF.date2num(datetime(year, 1, 1, 0, 0, 0), unitsRiver)
    # tendRiver = netCDF.date2num(datetime(year+1, 1, 1, 0, 0, 0), unitsRiver)
    # start and end indices in time for river discharge
    itstartRiver = bisect.bisect_left(datesRiver, datetime(year, 1, 1, 0, 0, 0))
    itendRiver = bisect.bisect_left(datesRiver, datetime(year+1, 1, 1, 0, 0, 0))
        axr.plot(tRiver[itstartRiver:itriver], Q[itstartRiver:itriver], '-', color='0.4')
        axr.plot(tRiver[itriver:itendRiver+1], Q[itriver:itendRiver+1], '-', color='0.4', alpha=0.3)


    ## Wind forcing ##
    w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
    # Wind time period to use
    unitsWind = w.variables['time'].units
    datesWind = netCDF.num2date(w.variables['time'][:], unitsWind)
        Uwind = w.variables['Uwind'][itwind,:,:]
        Vwind = w.variables['Vwind'][itwind,:,:]




# write table to file
with open('test_file.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in table]
