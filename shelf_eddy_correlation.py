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
import dateutil.relativedelta


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


years = np.arange(2004,2011)
months = np.arange(6,10)

## River forcing ##
r = netCDF.Dataset('/atch/raid1/zhangxq/Projects/txla_nesting6/TXLA_river_4dyes_2011.nc')
# River timing
unitsRiver = r.variables['river_time'].units
datesRiver = netCDF.num2date(r.variables['river_time'][:], unitsRiver)
tRiver = r.variables['river_time'][:]
# all of river input
Q = np.abs(r.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
r.close()

grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
gridnc = netCDF.Dataset(grid_filename)
iwind = gridnc.variables['lat_rho'][:]>27.5

# transport calculated previously
T = np.load('calcs/shelfconn/interannual-summertransportsum.npz')['transport']

table = np.empty((years.size,7*months.size+1))

# Loop through years and write a data point in each column for the row each the year
for j,year in enumerate(years):

    # Instantaneous discharge at the start of each month, Qi
    Qi = np.empty(months.size)
    for i,month in enumerate(months):
        itime = find(datesRiver<datetime(year,month,1,0,0,0))[-1] # starting index
        Qi[i] = Q[itime]
        # ti[i] = tRiver[itime] # time numbers for reference

    # Discharge cumulative from beginning of year through the beginning of each month, Qcum
    Qcum = np.empty(months.size)
    for i,month in enumerate(months):
        istart = find(datesRiver<datetime(year,1,1,0,0,0))[-1] # starting index, January
        iend = find(datesRiver<datetime(year,month,1,0,0,0))[-1] # ending index, current month
        Qcum[i] = Q[istart:iend].sum()

    # Discharge running sum over three months, Q3
    Q3 = np.empty(months.size)
    for i,month in enumerate(months):
        # datetime three months ago
        dt3 = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=3)
        istart = find(datesRiver<dt3)[-1] # starting index
        iend = find(datesRiver<datetime(year,month,1,0,0,0))[-1] # ending index, current month
        Q3[i] = Q[istart:iend].sum()



    ## Wind forcing: averaged over the broad shelf region and taken over two months ##
    baseloc = '/atch/raid1/zhangxq/Projects/narr_txla/'
    # if year==2010:
    #     w = netCDF.MFDataset([baseloc + 'txla_blk_narr_2009.nc', baseloc + 'txla_blk_narr_2010.nc'])
    # else:
    # w = netCDF.MFDataset(np.sort(glob.glob(baseloc + 'txla_blk_narr_200[' + str(year-1)[-1] + '-' + str(year)[-1] + '].nc')))
    w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
    # Wind time period to use
    unitsWind = w.variables['time'].units
    tWind = w.variables['time'][:]
    datesWind = netCDF.num2date(tWind, unitsWind)
    Wsmean = np.empty(months.size)
    Wsvar = np.empty(months.size)
    Wdmean = np.empty(months.size)
    Wdvar = np.empty(months.size)
    for i,month in enumerate(months):
        # datetime three months ago
        dt2 = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=2)
        istart = find(datesWind<dt2)[-1] # starting index, 2 months ago
        iend = find(datesWind<datetime(year,month,1,0,0,0))[-1] # ending index, current month

        Uwind = w.variables['Uwind'][istart:iend,:,:]
        Uwind = Uwind[:,iwind]
        Vwind = w.variables['Vwind'][istart:iend,:,:]
        Vwind = Vwind[:,iwind]
        
        Swind = np.sqrt(Uwind**2 + Vwind**2)
        Uwindmean = Uwind.mean()
        Vwindmean = Vwind.mean()        

        # Wind magnitude mean, Wsmean
        Wsmean[i] = Swind.mean()

        # Wind magnitude var, Wsvar
        Wsvar[i] = np.var(Swind)

        # Wind direction mean, Wdmean
        Wdmean[i] = np.rad2deg(np.arctan2(Vwindmean, Uwindmean))

        # Wind direction var, Wdvar
        Wdvar[i] = np.var(np.arctan2(Vwind, Uwind))

    w.close()

    # Combine together into one row (actually put into matrix)
    table[j,:] = np.hstack((T[j], Qi, Qcum, Q3, Wsmean, Wsvar, Wdmean, Wdvar))


# Headers
headers = ('Transport', 'Instantaneous-discharge-Qi-Jun', 'Qi-Jul', 'Qi-Aug', 'Qi-Sep',
            'Cumulative-discharge-Qcum-Jun', 'Qcum-Jul', 'Qcum-Aug', 'Qcum-Sep',
            'Running-sum-discharge-Q3-Jun', 'Q3-Jul', 'Q3-Aug', 'Q3-Sep',
            'Wind-magnitude-mean-Wsmean-Jun', 'Wsmean-Jul', 'Wsmean-Aug', 'Wsmean-Sep',
            'Wind-magnitude-variance-Wsvar-Jun', 'Wsvar-Jul', 'Wsvar-Aug', 'Wsvar-Sep',
            'Wind-direction-mean-Wdmean-Jun', 'Wdmean-Jul', 'Wdmean-Aug', 'Wdmean-Sep',
            'Wind-direction-variance-Wdvar-Jun', 'Wdvar-Jul', 'Wdvar-Aug', 'Wdvar-Sep')
# headers = ('Transport', 'Instantaneous-discharge-Qi-Jan', 'Qi-Feb', 'Qi-Mar', 'Qi-Apr', 'Qi-May', 'Qi-Jun', 'Qi-Jul',
#             'Qi-Aug', 'Qi-Sep',
#             'Cumulative-discharge-Qcum-Feb', 'Qcum-Mar', 'Qcum-Apr', 'Qcum-May', 'Qcum-Jun', 'Qcum-Jul',
#             'Qcum-Aug', 'Qcum-Sep',
#             'Running-sum-discharge-Q3-Jan', 'Q3-Feb', 'Q3-Mar', 'Q3-Apr', 'Q3-May', 'Q3-Jun', 'Q3-Jul',
#             'Q3-Aug', 'Q3-Sep',
#             'Wind-magnitude-mean-Wsmean-Jan', 'Wsmean-Feb', 'Wsmean-Mar', 'Wsmean-Apr', 'Wsmean-May',
#             'Wsmean-Jun', 'Wsmean-Jul', 'Wsmean-Aug', 'Wsmean-Sep',
#             'Wind-magnitude-variance-Wsvar-Jan', 'Wsvar-Feb', 'Wsvar-Mar', 'Wsvar-Apr', 'Wsvar-May',
#             'Wsvar-Jun', 'Wsvar-Jul', 'Wsvar-Aug', 'Wsvar-Sep',
#             'Wind-direction-mean-Wdmean-Jan', 'Wdmean-Feb', 'Wdmean-Mar', 'Wdmean-Apr', 'Wdmean-May',
#             'Wdmean-Jun', 'Wdmean-Jul', 'Wdmean-Aug', 'Wdmean-Sep',
#             'Wind-direction-variance-Wdvar-Jan', 'Wdvar-Feb', 'Wdvar-Mar', 'Wdvar-Apr', 'Wdvar-May',
#             'Wdvar-Jun', 'Wdvar-Jul', 'Wdvar-Aug', 'Wdvar-Sep')

# write table to file
with open('table-new.txt', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(headers)
    [writer.writerow(r) for r in table]
