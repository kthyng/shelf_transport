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
# import tracpy
from datetime import datetime, timedelta
import glob
# from cmPong import cmPong
from matplotlib.mlab import find
import bisect
from matplotlib import delaunay
import csv
import dateutil.relativedelta
# import op
import scipy.stats
from calendar import monthrange
import seaborn as sns
sns.set(style='whitegrid', color_codes=True)
import pandas as pd
import cmocean

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


maketable = False
runcorr = True  # run correlation coefficients
doplot = True
explore = False  # True to do exploratory plots, False to do polished plot
whichseason = 'summer'  # 'winter' or 'summer' or 'dispersion'
makesmalltable = False  # create a subset table for analysis in R
if whichseason == 'dispersion':
    col = 0
elif whichseason == 'winter':
    col = 1
elif whichseason == 'summer':
    col = 2
table = 'table-sustr-100m-region34'
without2008 = False
makeplot = False  # for shelf_transport paper

# If doing small table creation, need a list of things to include
if makesmalltable:
    if whichseason == 'winter':
        Vars = ['Transport_winter_R3', 'Wdm3-Apr', 'Wdm3-Nov_previous', 'Q1-Dec_previous', 'Wsm3-Nov_previous']
        name = table + '-winter'
    elif whichseason == 'summer':
        Vars = ['Transport_summer_R4', 'Q1-Apr', 'Qi-Aug_previous', 'Wsm1-Jul', 'Wsm1-Oct', 'Wsm3-Aug_previous']
        name = table + '-summer'


headers = ('Dispersion_summer_R4', 'Transport_winter_R3', 'Transport_summer_R4',
            'Instantaneous-discharge-Qi-Jan', 'Qi-Feb', 'Qi-Mar', 'Qi-Apr', 'Qi-May', 'Qi-Jun', 'Qi-Jul',
            'Qi-Aug', 'Qi-Sep', 'Qi-Oct', 'Qi-Nov', 'Qi-Dec',
            'Cumulative-discharge-Qcum-Jan', 'Qcum-Feb', 'Qcum-Mar', 'Qcum-Apr', 'Qcum-May', 'Qcum-Jun', 'Qcum-Jul',
            'Qcum-Aug', 'Qcum-Sep', 'Qcum-Oct', 'Qcum-Nov', 'Qcum-Dec',
            'Running-sum-discharge-Q1-Jan', 'Q1-Feb', 'Q1-Mar', 'Q1-Apr', 'Q1-May', 'Q1-Jun', 'Q1-Jul',
            'Q1-Aug', 'Q1-Sep', 'Q1-Oct', 'Q1-Nov', 'Q1-Dec',
            'Running-sum-discharge-Q2-Jan', 'Q2-Feb', 'Q2-Mar', 'Q2-Apr', 'Q2-May', 'Q2-Jun', 'Q2-Jul',
            'Q2-Aug', 'Q2-Sep', 'Q2-Oct', 'Q2-Nov', 'Q2-Dec',
            'Running-sum-discharge-Q3-Jan', 'Q3-Feb', 'Q3-Mar', 'Q3-Apr', 'Q3-May', 'Q3-Jun', 'Q3-Jul',
            'Q3-Aug', 'Q3-Sep', 'Q3-Oct', 'Q3-Nov', 'Q3-Dec',
            'Wind-magnitude-mean-Wsm1-Jan', 'Wsm1-Feb', 'Wsm1-Mar', 'Wsm1-Apr', 'Wsm1-May',
            'Wsm1-Jun', 'Wsm1-Jul', 'Wsm1-Aug', 'Wsm1-Sep', 'Wsm1-Oct', 'Wsm1-Nov', 'Wsm1-Dec',
            'Wind-magnitude-mean-Wsm2-Jan', 'Wsm2-Feb', 'Wsm2-Mar', 'Wsm2-Apr', 'Wsm2-May',
            'Wsm2-Jun', 'Wsm2-Jul', 'Wsm2-Aug', 'Wsm2-Sep', 'Wsm2-Oct', 'Wsm2-Nov', 'Wsm2-Dec',
            'Wind-magnitude-mean-Wsm3-Jan', 'Wsm3-Feb', 'Wsm3-Mar', 'Wsm3-Apr', 'Wsm3-May',
            'Wsm3-Jun', 'Wsm3-Jul', 'Wsm3-Aug', 'Wsm3-Sep', 'Wsm3-Oct', 'Wsm3-Nov', 'Wsm3-Dec',
            'Wind-magnitude-variance-Wsv1-Jan', 'Wsv1-Feb', 'Wsv1-Mar', 'Wsv1-Apr', 'Wsv1-May',
            'Wsv1-Jun', 'Wsv1-Jul', 'Wsv1-Aug', 'Wsv1-Sep', 'Wsv1-Oct', 'Wsv1-Nov', 'Wsv1-Dec',
            'Wind-magnitude-variance-Wsv2-Jan', 'Wsv2-Feb', 'Wsv2-Mar', 'Wsv2-Apr', 'Wsv2-May',
            'Wsv2-Jun', 'Wsv2-Jul', 'Wsv2-Aug', 'Wsv2-Sep', 'Wsv2-Oct', 'Wsv2-Nov', 'Wsv2-Dec',
            'Wind-magnitude-variance-Wsv3-Jan', 'Wsv3-Feb', 'Wsv3-Mar', 'Wsv3-Apr', 'Wsv3-May',
            'Wsv3-Jun', 'Wsv3-Jul', 'Wsv3-Aug', 'Wsv3-Sep', 'Wsv3-Oct', 'Wsv3-Nov', 'Wsv3-Dec',
            'Wind-direction-mean-Wdm1-Jan', 'Wdm1-Feb', 'Wdm1-Mar', 'Wdm1-Apr', 'Wdm1-May',
            'Wdm1-Jun', 'Wdm1-Jul', 'Wdm1-Aug', 'Wdm1-Sep', 'Wdm1-Oct', 'Wdm1-Nov', 'Wdm1-Dec',
            'Wind-direction-mean-Wdm2-Jan', 'Wdm2-Feb', 'Wdm2-Mar', 'Wdm2-Apr', 'Wdm2-May',
            'Wdm2-Jun', 'Wdm2-Jul', 'Wdm2-Aug', 'Wdm2-Sep', 'Wdm2-Oct', 'Wdm2-Nov', 'Wdm2-Dec',
            'Wind-direction-mean-Wdm3-Jan', 'Wdm3-Feb', 'Wdm3-Mar', 'Wdm3-Apr', 'Wdm3-May',
            'Wdm3-Jun', 'Wdm3-Jul', 'Wdm3-Aug', 'Wdm3-Sep', 'Wdm3-Oct', 'Wdm3-Nov', 'Wdm3-Dec',
            'Wind-direction-variance-Wdv1-Jan', 'Wdv1-Feb', 'Wdv1-Mar', 'Wdv1-Apr', 'Wdv1-May',
            'Wdv1-Jun', 'Wdv1-Jul', 'Wdv1-Aug', 'Wdv1-Sep', 'Wdv1-Oct', 'Wdv1-Nov', 'Wdv1-Dec',
            'Wind-direction-variance-Wdv2-Jan', 'Wdv2-Feb', 'Wdv2-Mar', 'Wdv2-Apr', 'Wdv2-May',
            'Wdv2-Jun', 'Wdv2-Jul', 'Wdv2-Aug', 'Wdv2-Sep', 'Wdv2-Oct', 'Wdv2-Nov', 'Wdv2-Dec',
            'Wind-direction-variance-Wdv3-Jan', 'Wdv3-Feb', 'Wdv3-Mar', 'Wdv3-Apr', 'Wdv3-May',
            'Wdv3-Jun', 'Wdv3-Jul', 'Wdv3-Aug', 'Wdv3-Sep', 'Wdv3-Oct', 'Wdv3-Nov', 'Wdv3-Dec')


def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

# grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid_filename = '../grid.nc'
dtemp = netCDF.Dataset(grid_filename)
anglev = dtemp.variables['angle'][:]
dtemp.close()

if makesmalltable:
    d = np.loadtxt(table + '.txt', skiprows=1)

    # double the headers for the correlations
    headers2 = list(headers)
    for i in xrange(len(headers2)):
        headers2[i] += '_previous'
    headers = np.hstack((headers, headers2))

    table = np.empty((11, len(Vars))) # Vars metrics being included and also 1 transport

    for j, Var in enumerate(Vars):
        ind = list(headers).index(Var) # find index of var
        if 'previous' in Var:  # need to shift accordingly
            table[:,j] = d[:-1, ind-206]
        else:  # need to shift accordingly
            table[:,j] = d[1:, ind]

    # write table to file
    with open(name + '.txt', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(Vars)
        [writer.writerow(r) for r in table]



elif runcorr:
    d = np.loadtxt(table + '.txt', skiprows=1)

    without2008 = False

    if without2008:
        # skip 2008
        d2 = np.vstack((d[0:4,:], d[5:,:]))
        d = d2.copy()

    # which column to use for transport; skip fake 2003
    T = d[1:, col]

    # if whichseason == 'winter':
    #     T = d[:,0] # the transport is the first column
    # elif whichseason == 'summer':
    #     T = d[:,1] # the transport is the first column
    r = np.empty(d.shape[1]*2)  # to store the correlation coefficients
    p = np.empty(d.shape[1]*2)  # to store the correlation coefficients
    for i in range(d.shape[1]):
        # metric to compare with, 2004-2014, year before transport
        # store in the immediate location
        comp = d[1:, i]
        ind = ~np.isnan(comp)
        r[i], p[i] = scipy.stats.pearsonr(T[ind], comp[ind])
        # metric to compare with, 2003-2013, year before transport
        # store in the index+d.shape[1]
        comp = d[:-1, i]
        ind = ~np.isnan(comp)
        r[i+d.shape[1]], p[i+d.shape[1]] = scipy.stats.pearsonr(T[ind], comp[ind])

    # double the headers for the correlations
    headers2 = list(headers)
    for i in range(len(headers2)):
        headers2[i] += '_previous'
    headers = np.hstack((headers, headers2))

    ind = np.argsort(p)  # indices that would sort p in ascending order
    # print the best ones
    for i in range(20):
        print(headers[ind[i]] + ': p=%1.4f, r=%1.2f' % (p[ind[i]], r[ind[i]]))

    # hard wiring this for now. Sorry future Kristen!
    if whichseason == 'winter':
        ind = list(headers).index('Wdm3-Apr')
        Wbest = d[1:, ind]
        Wbestr = r[ind]
        Wbestp = p[ind]
    elif whichseason == 'summer':
        ind = list(headers).index('Q1-Apr')
        Qbest = d[1:, ind]/1e4
        Qbestr = r[ind]
        Qbestp = p[ind]

    # # save best entries to plot
    # # Top river one
    # inan = ~(np.isnan(r)) * ~(r==1)
    # ind = np.argmax(abs(r[inan]))
    # while ('Qbest' not in locals()):
    #     if 'Q' in np.asarray(headers)[inan][ind]:
    #         Qbestr = r[inan][ind]
    #         Qbestp = p[inan][ind]
    #         import pdb; pdb.set_trace()
    #         Qbest = d[:,inan][:,ind]
    #         Qbestname = np.asarray(headers)[inan][ind]
    #     inan *= ~(r==r[inan][ind])
    #     ind = np.argmax(r[inan])

    # # Top wind one
    # inan = ~(np.isnan(r)) * ~(r==1)
    # ind = np.argmax(abs(r[inan]))
    # while ('Wbest' not in locals()):
    #     if 'W' in np.asarray(headers)[inan][ind]:
    #         Wbestr = r[inan][ind]
    #         Wbestp = p[inan][ind]
    #         Wbest = d[:,inan][:,ind]
    #         Wbestname = np.asarray(headers)[inan][ind]
    #     inan *= ~(r==r[inan][ind])
    #     ind = np.argmax(r[inan])


# use info from runcorr
if doplot:

    if without2008:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mappable = ax.scatter(Qbest, Wbest, c=T, cmap='Reds', s=300)
        ax.set_xlabel('July river discharge')
        ax.set_ylabel('July-August wind direction variance')
        cb = fig.colorbar(mappable)
        cb.set_label('Relative cross-shelf transport')
        fig.savefig('figures/summer-transport-correlations-without2008.pdf', bbox_inches='tight')
        fig.savefig('figures/summer-transport-correlations-without2008.png', bbox_inches='tight')
    else:
        if explore:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            mappable = ax.scatter(Qbest, Wbest, c=T, cmap='Reds', s=300)
            ax.set_xlabel(Qbestname + ', r=' + str(Qbestr) + ', p=' + str(Qbestp))
            ax.set_ylabel(Wbestname + ', r=' + str(Wbestr) + ', p=' + str(Wbestp))
            cb = fig.colorbar(mappable)
            cb.set_label('Relative cross-shelf transport')
            # fig.savefig('figures/' + whichseason + '-transport-correlations.pdf', bbox_inches='tight')
            # fig.savefig('figures/' + whichseason + '-transport-correlations.png', bbox_inches='tight')

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(T, Qbest, 'o')
            ax.set_xlabel('Transport')
            ax.set_ylabel(Qbestname + ', r=' + str(Qbestr) + ', p=' + str(Qbestp))

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(T, Wbest, 'o')
            ax.set_xlabel('Transport')
            ax.set_ylabel(Wbestname + ', r=' + str(Wbestr) + ', p=' + str(Wbestp))

        else:
            if whichseason == 'winter':
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111)
                # ax.plot(T, Wbest, 'o', ms=10, color='0.2', mec='k')
                ax.set_xlim(T.min()-0.1, T.max()+0.1)
                ax.set_ylim(Wbest.min()-5, Wbest.max()+5)
                sns.regplot(x=T, y=Wbest, ax=ax, scatter_kws={'s': 50});
                ax.set_xlabel('Transport relative to mean [%]', fontsize=14)
                ax.set_ylabel('January-February-March mean wind direction [deg]', fontsize=14)
                plt.tick_params('both', labelsize=14)
                ax.text(0.1, 0.15, 'r=%1.2f' % Wbestr, color='#15317E', transform=ax.transAxes, alpha=0.7, fontsize=24)
                ax.text(0.1, 0.07, 'p=%1.4f' % Wbestp, color='#15317E', transform=ax.transAxes, alpha=0.7, fontsize=24)
                # ax.set_frame_on(False)
                fig.savefig('figures/winter-transport-correlations.pdf', bbox_inches='tight')
            elif whichseason == 'summer':
                # # with two variables
                # fig = plt.figure(figsize=(6,6))
                # ax = fig.add_subplot(111)
                # mappable = ax.scatter(d[:,29]/1e5, d[:,117], c=T, cmap='Reds', s=300)
                # # ax.set_xlim(T.min()-0.1, T.max()+0.1)
                # ax.set_xlabel('March river discharge/10^5 [m$^3\!$/s]')
                # ax.set_ylabel('June wind speed variance [m/s]')
                # ax.text(0.75, 0.15, '$r^2_{both}$=0.81', color='r', transform=ax.transAxes, alpha=0.7)
                # ax.text(0.75, 0.08, '$r^2_{river}$=0.78', color='r', transform=ax.transAxes, alpha=0.7)
                # cb = fig.colorbar(mappable)
                # cb.set_label('Transport relative to mean [%]')
                # ax.set_frame_on(False)
                # fig.savefig('figures/summer-transport-correlations2.pdf', bbox_inches='tight')

                # with one variable
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111)
                # mappable = ax.plot(T, d[:,29]/1e4, 'o', ms=10, color='0.2', mec='k')
                ax.set_xlim(T.min()-0.1, T.max()+0.1)
                ax.set_ylim(1.0, 2.8)
                sns.regplot(x=T, y=Qbest, ax=ax, scatter_kws={'s': 50});
                ax.set_xlabel('Transport relative to mean [%]', fontsize=14)
                ax.set_ylabel('Mean March river discharge~' + r'$\left[10^4\mathrm{m}^3\mathrm{s}^{-1}\right]$', fontsize=14)
                plt.tick_params('both', labelsize=14)
                ax.text(0.6, 0.13, 'r=%1.2f' % Qbestr, color='#15317E', transform=ax.transAxes, alpha=0.7, fontsize=24)
                ax.text(0.6, 0.04, 'p=%1.4f' % Qbestp, color='#15317E', transform=ax.transAxes, alpha=0.7, fontsize=24)
                # ax.set_frame_on(False)
                fig.savefig('figures/summer-transport-correlations1.pdf', bbox_inches='tight')


if maketable:

    # for years, include year before simulations so that I can look at metrics from beforehand
    # for all seasons
    years = np.arange(2014, 2015)
    months = np.arange(1, 13)

    # Model output ##
    # currents_filenames = []
    # currents_filenames.extend(np.sort(glob.glob('/home/kthyng/shelf/200[3-9]/ocean_his_????.nc')))
    # currents_filenames.extend(np.sort(glob.glob('/home/kthyng/shelf/201[0-3]/ocean_his_????.nc')))
    # currents_filenames.extend((glob.glob('/home/kthyng/shelf/2014/ocean_his_??.nc')))
    # m = netCDF.MFDataset(currents_filenames)

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    nc = netCDF.Dataset(loc)  # for wind
    t = nc.variables['ocean_time']
    datesw = netCDF.num2date(t[:], t.units)

    # Surface stress
    sustr = nc.variables['sustr']; svstr = nc.variables['svstr']

    ## River forcing ##
    # r1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/txla_nesting6/TXLA_river_4dyes_2011.nc')
    r1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012.nc')  # use for through 2011
    r2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012_2014.nc')  # use for 2012-2014
    # River timing
    tr1 = r1.variables['river_time']
    tunitsr1 = tr1.units
    # interpolate times for this data file since at the 12 hours mark instead of beginning of the day
    tr1 = op.resize(tr1, 0)
    datesr1 = netCDF.num2date(tr1[:], tunitsr1)
    tr2 = r2.variables['river_time']
    datesr2 = netCDF.num2date(tr2[:], tr2.units)
    # all of river input
    Q1 = np.abs(r1.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
    # interpolate this like for time
    Q1 = op.resize(Q1, 0)
    Q2 = np.abs(r2.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
    # pdb.set_trace()

    # Combine river info into one dataset
    iend1 = find(datesr1<datetime(2012,1,1,0,0,0))[-1] # ending index for file 1
    tr = np.concatenate((tr1[:iend1], tr2[:]), axis=0)
    datesr = np.concatenate((datesr1[:iend1], datesr2))
    Q = np.concatenate((Q1[:iend1], Q2))
    r1.close(); r2.close()

    # grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    # gridnc = netCDF.Dataset(grid_filename)
    iwind = nc.variables['lat_rho'][:] > 27.5
    iwindu = nc.variables['lat_u'][:] > 27.5
    iwindv = nc.variables['lat_v'][:] > 27.5

    # transport calculated previously
    Tw = np.load('calcs/shelfconn/interannual-winter3transportsum.npz')['transport']
    Ts = np.load('calcs/shelfconn/interannual-summer4transportsum.npz')['transport']
    Tw = np.hstack((0, Tw))  # add a spot for 2003, which is run up time
    Ts = np.hstack((0, Ts))  # add a spot for 2003, which is run up time

    # integrate dispersion in box
    # ADD TO TABLE FILE
    D2s = np.load('calcs/dispersion/interannual-summer4dispersionsum.npz')['transport']
    D2s = np.hstack((0, D2s))  # add a spot for 2003, which is run up time

    table = np.empty((years.size,17*months.size+2)) # 17 metrics being calculated and also 2 transports

    # Loop through years and write a data point in each column for the row each the year
    for j,year in enumerate(years):
        print(year)

        # if year == 2014:
        #     months = np.arange(1,10)

        # Instantaneous discharge at the start of each month, Qi
        Qi = np.zeros(months.size)
        for i,month in enumerate(months):
            itime = find(datesr<datetime(year,month,1,0,0,0))[-1] # starting index
            Qi[i] = Q[itime]
            # ti[i] = tr1[itime] # time numbers for reference

        # Discharge cumulative from beginning of year through the beginning of each month, Qcum
        Qcum = np.zeros(months.size)
        for i,month in enumerate(months):
            istart = find(datesr<datetime(year,1,1,0,0,0))[-1] # starting index, January
            iend = find(datesr<datetime(year,month,1,0,0,0))[-1] # ending index, current month
            Qcum[i] = Q[istart:iend].sum()

        # Discharge running mean over 1 month, Q1
        Q1 = np.zeros(months.size)
        for i, month in enumerate(months):
            # datetime 1 month ago
            dt = datetime(year, month, 1, 0, 0, 0) - dateutil.relativedelta.relativedelta(months=1)
            istart = find(datesr < dt)[-1]  # starting index
            iend = find(datesr < datetime(year, month, 1, 0, 0, 0))[-1]  # ending index, current month
            ndays = monthrange(dt.year, dt.month)[1]  # number of days in the month before
            Q1[i] = Q[istart:iend].sum()/ndays

        # Discharge running mean over 2 months, Q2
        Q2 = np.zeros(months.size)
        for i, month in enumerate(months):
            # datetime 2 months ago
            dt = datetime(year, month, 1, 0, 0, 0) - dateutil.relativedelta.relativedelta(months=2)
            istart = find(datesr < dt)[-1]  # starting index
            iend = find(datesr < datetime(year, month, 1, 0, 0, 0))[-1]  # ending index, current month
            ndays = monthrange(dt.year, dt.month)[1]  # days from 1st month
            dt = dt + dateutil.relativedelta.relativedelta(months=1)  # just to step forward for ndays
            ndays += monthrange(dt.year, dt.month)[1]  # days from 2nd month
            Q2[i] = Q[istart:iend].sum()/ndays

        # Discharge running mean over three months, Q3
        Q3 = np.zeros(months.size)
        for i, month in enumerate(months):
            # datetime three months ago
            dt = datetime(year, month, 1, 0, 0, 0) - dateutil.relativedelta.relativedelta(months=3)
            istart = find(datesr < dt)[-1]  # starting index
            iend = find(datesr < datetime(year, month, 1, 0, 0, 0))[-1]  # ending index, current month
            ndays = monthrange(dt.year, dt.month)[1]  # days from 1st month
            dt = dt + dateutil.relativedelta.relativedelta(months=1)  # just to step forward for ndays
            ndays += monthrange(dt.year, dt.month)[1]  # days from 2nd month
            dt = dt + dateutil.relativedelta.relativedelta(months=1)  # just to step forward for ndays
            ndays += monthrange(dt.year, dt.month)[1]  # days from 3rd month
            Q3[i] = Q[istart:iend].sum()/ndays

        ## Wind forcing: averaged over the broad shelf region and taken over 1, 2, and 3 months ##
        # # There are multiple file locations
        # if year <= 2012:
        #     # have to fudge for the metrics early in 2003 that should be calculated from late 2002,
        #     # but should be ok since nothing should really depend on those things.
        #     if year == 2003:
        #         w1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
        #     else:
        #         w1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year-1) + '.nc')
        #     w2 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
        # elif year == 2013:
        #     w1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_2012.nc')
        #     w2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
        # elif year == 2014:
        #     w1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
        #     w2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2014.nc')

        # # Wind timing
        # # pdb.set_trace()
        # tw1 = w1.variables['time']
        # if year == 2003:  # fake the year
        #     datesw1 = netCDF.num2date(tw1[:], tw1.units.replace('/', '-')) - \
        #                  dateutil.relativedelta.relativedelta(years=1)
        # else:
        #     datesw1 = netCDF.num2date(tw1[:], tw1.units.replace('/', '-'))
        # tw2 = w2.variables['time']
        # datesw2 = netCDF.num2date(tw2[:], tw2.units.replace('/', '-'))

        # # Combine wind info here since files can't be combined in MFDataset
        # Uwind1 = w1.variables['Uwind'][:]
        # Vwind1 = w1.variables['Vwind'][:]
        # # rotate onto cartesian grid
        # Uwind1, Vwind1 = rot2d(Uwind1, Vwind1, anglev)
        # Uwind2 = w2.variables['Uwind'][:]
        # Vwind2 = w2.variables['Vwind'][:]
        # # rotate onto cartesian grid
        # Uwind2, Vwind2 = rot2d(Uwind2, Vwind2, anglev)
        # Uwind = np.concatenate((Uwind1, Uwind2[1:, :, :]), axis=0)
        # import pdb; pdb.set_trace()
        # del(Uwind1, Uwind2)
        # Uwind = Uwind[:, iwind]  # only use winds from certain area of domain

        # Vwind = np.concatenate((Vwind1, Vwind2[1:, :, :]), axis=0)
        # del(Vwind1, Vwind2)
        # Vwind = Vwind[:, iwind]  # only use winds from certain area of domain
        # Swind = np.sqrt(Uwind**2 + Vwind**2)

        # tw = np.concatenate((tw1, tw2[1:]), axis=0)
        # datesw = np.concatenate((datesw1, datesw2[1:]), axis=0)
        # w1.close(); w2.close()

        Wsm1 = np.zeros(months.size); Wsv1 = np.zeros(months.size)
        Wdm1 = np.zeros(months.size); Wdv1 = np.zeros(months.size)
        Wsm2 = np.zeros(months.size); Wsv2 = np.zeros(months.size)
        Wdm2 = np.zeros(months.size); Wdv2 = np.zeros(months.size)
        Wsm3 = np.zeros(months.size); Wsv3 = np.zeros(months.size)
        Wdm3 = np.zeros(months.size); Wdv3 = np.zeros(months.size)
        for i, month in enumerate(months):

            if (year == 2003) and ((month == 1) or (month == 2) or (month == 3) or (month == 4) or (month == 5)):
                continue  # don't have info for here
            elif (year == 2014) and (month >= 9):
                continue  # don't have info for here

            # datetime 1 month ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=1)
            istart1 = find(datesw<dt)[-1] # starting index
            iend1 = find(datesw<datetime(year,month,1,0,0,0))[-1] # ending index, current month
            # datetime 2 months ago
            # pdb.set_trace()
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=2)
            istart2 = find(datesw<dt)[-1] # starting index
            # do each month incrementally since not enough memory otherwise
            iend2 = find(datesw<(datetime(year,month,1,0,0,0)-dateutil.relativedelta.relativedelta(months=1)))[-1] # ending index, current month
            # datetime 3 months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=3)
            istart3 = find(datesw<dt)[-1] # starting index
            iend3 = find(datesw<(datetime(year,month,1,0,0,0)-dateutil.relativedelta.relativedelta(months=2)))[-1] # ending index, current month

            # 1st month
            # if month==2:
            #     pdb.set_trace()
            sustr1 = sustr[istart1:iend1, :]; svstr1 = svstr[istart1:iend1, :]
            sustr1 = op.resize(sustr1, 2)[:,1:-1,:]
            svstr1 = op.resize(svstr1, 1)[:,:,1:-1]
            # rotate onto cartesian grid
            sustr1, svstr1 = rot2d(sustr1, svstr1, anglev[1:-1,1:-1])
            sustr1 = sustr1[:, iwind[1:-1,1:-1]]
            svstr1 = svstr1[:, iwind[1:-1,1:-1]]
            # sum for later
            l1 = sustr1.shape[0]  # to divide by later for mean
            ssstrsum = np.sqrt(sustr1**2 + svstr1**2).sum()
            sustrsum = sustr1.sum()
            svstrsum = svstr1.sum()
            del(sustr1,svstr1)
            # pdb.set_trace()
            Wsm1[i] = ssstrsum/l1  # Wind magnitude mean, Wsm
            Wdm1[i] = np.rad2deg(np.arctan2(svstrsum/l1, sustrsum/l1))  # Wind direction mean, Wdm

            # 2nd month
            sustr2 = sustr[istart2:iend2, :]; svstr2 = svstr[istart2:iend2, :]
            sustr2 = op.resize(sustr2, 2)[:,1:-1,:]
            svstr2 = op.resize(svstr2, 1)[:,:,1:-1]
            # rotate onto cartesian grid
            sustr2, svstr2 = rot2d(sustr2, svstr2, anglev[1:-1,1:-1])
            sustr2 = sustr2[:, iwind[1:-1,1:-1]]
            svstr2 = svstr2[:, iwind[1:-1,1:-1]]
            # sum for later
            l2 = sustr2.shape[0]  # to divide by later for mean
            ssstrsum += np.sqrt(sustr2**2 + svstr2**2).sum()
            sustrsum += sustr2.sum()
            svstrsum += svstr2.sum()
            del(sustr2,svstr2)
            Wsm2[i] = ssstrsum/(l1 + l2)  # Wind magnitude mean, Wsm
            Wdm2[i] = np.rad2deg(np.arctan2(svstrsum/(l1+l2), sustrsum/(l1+l2)))  # Wind direction mean, Wdm

            # 3rd month
            sustr3 = sustr[istart3:iend3, :]; svstr3 = svstr[istart3:iend3, :]
            sustr3 = op.resize(sustr3, 2)[:,1:-1,:]
            svstr3 = op.resize(svstr3, 1)[:,:,1:-1]
            # rotate onto cartesian grid
            sustr3, svstr3 = rot2d(sustr3, svstr3, anglev[1:-1,1:-1])
            sustr3 = sustr3[:, iwind[1:-1,1:-1]]
            svstr3 = svstr3[:, iwind[1:-1,1:-1]]
            # sum for later
            l3 = sustr3.shape[0]  # to divide by later for mean
            ssstrsum += np.sqrt(sustr3**2 + svstr3**2).sum()
            sustrsum += sustr3.sum()
            svstrsum += svstr3.sum()
            del(sustr3,svstr3)
            Wsm3[i] = ssstrsum/(l1 + l2 + l3)  # Wind magnitude mean, Wsm
            Wdm3[i] = np.rad2deg(np.arctan2(svstrsum/(l1+l2+l3), sustrsum/(l1+l2+l3)))  # Wind direction mean, Wdm
            # pdb.set_trace()
            # # Wind magnitude var, Wsv
            # Wsv1[i] = np.var(ssstr1)
            # Wsv2[i] = np.var(ssstr2)
            # Wsv3[i] = np.var(ssstr3)

            # # Wind direction var, Wdv
            # ang = np.arctan2(svstr1, sustr1)
            # # unwrap angles before taking the variance
            # ang = np.rad2deg(np.unwrap(np.deg2rad(ang)))
            # Wdv1[i] = np.var(ang)

            # ang = np.arctan2(svstr2, sustr2)
            # ang = np.rad2deg(np.unwrap(np.deg2rad(ang)))
            # Wdv2[i] = np.var(ang)

            # ang = np.arctan2(svstr3, sustr3)
            # ang = np.rad2deg(np.unwrap(np.deg2rad(ang)))
            # Wdv3[i] = np.var(ang)

        # Unwrap mean directions
        Wdm1 = np.rad2deg(np.unwrap(np.deg2rad(Wdm1)))
        Wdm2 = np.rad2deg(np.unwrap(np.deg2rad(Wdm2)))
        Wdm3 = np.rad2deg(np.unwrap(np.deg2rad(Wdm3)))

        # Change to compass coordinates instead of geometric and make positive
        # between 0 and 360
        Wdm1 = -Wdm1 + 90  # change to compass
        ind = Wdm1 < 0  # get between 0 and 360
        Wdm1[ind] = Wdm1[ind] + 360
        Wdm2 = -Wdm2 + 90  # change to compass
        ind = Wdm2 < 0  # get between 0 and 360
        Wdm2[ind] = Wdm2[ind] + 360
        Wdm3 = -Wdm3 + 90  # change to compass
        ind = Wdm3 < 0  # get between 0 and 360
        Wdm3[ind] = Wdm3[ind] + 360


        # Combine together into one row (actually put into matrix)
        table[j,:] = np.hstack((Tw[j], Ts[j], Qi, Qcum, Q1, Q2, Q3, Wsm1, Wsm2, Wsm3,
            Wsv1, Wsv2, Wsv3, Wdm1, Wdm2, Wdm3, Wdv1, Wdv2, Wdv3))

        # write table to file - to have something in case it breaks
        with open(table + str(year) + '.txt', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(headers)
            [writer.writerow(r) for r in table]

    # write table to file
    with open(table + '.txt', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(headers)
        [writer.writerow(r) for r in table]

# # Bootstrapping example from Rob
#  rr, pp = [], []
#  for n in range(100000):
#    .....:     ind = np.random.randint(0, 11, size=11)
#    .....:     r, p = scipy.stats.pearsonr(T[ind], Wbest[ind])
#    .....:     rr.append(r)
#    .....:     pp.append(p)
#    .....:
# hist(rr, bins=100)100)


if makeplot:

    # list of conditionals of headers that are wind-based
    Wboolsp = [('W' in header) and ('previous' in header) and (not 'Wsv' in header) and (not 'Wdv' in header) for header in headers]
    Windsp = find(Wboolsp)  # indices
    # Wheadersp = headers[Windsp]  # Wind-based headers
    Wbools = [('W' in header) and (not 'previous' in header) and (not 'Wsv' in header) and (not 'Wdv' in header) for header in headers]
    Winds = find(Wbools)  # indices
    Wheaders = headers[Winds]  # Wind-based headers
    # list of conditionals of headers that are river-based
    Qboolsp = [('Q' in header) and ('previous' in header) for header in headers]
    # Qboolsp = [('Q' in header) and ('previous' in header) and (not 'Cumulative-discharge' in header) for header in headers]
    Qindsp = find(Qboolsp)  # indices
    # Qheadersp = headers[Qindsp]  # River-based headers
    Qbools = [('Q' in header) and (not 'previous' in header) for header in headers]
    # Qbools = [('Q' in header) and (not 'previous' in header) and (not 'Cumulative-discharge' in header) for header in headers]
    Qinds = find(Qbools)  # indices
    Qheaders = headers[Qinds]  # River-based headers


    # figure()
    # scatter(abs(r[Winds]), p[Winds], s=500, facecolors='none', edgecolors='b', linewidths=2, alpha=0.5)
    # plot(abs(r[Rinds]), p[Rinds], 'o', color='r', ms=15, alpha=0.5)

    # previous year
    rwinterp = np.hstack((abs(r[Windsp]), abs(r[Qindsp])))[:, np.newaxis]  # all valid r values for winter, previous year
    pwinterp = np.hstack((abs(p[Windsp]), abs(p[Qindsp])))[:, np.newaxis]
    # pwinterp = -(1-pwinterp)  # to make them look different in colormap
    # concurrent year
    rwinter = np.hstack((abs(r[Winds]), abs(r[Qinds])))[:, np.newaxis]  # all valid r values for winter
    pwinter = np.hstack((abs(p[Winds]), abs(p[Qinds])))[:, np.newaxis]
    # pwinter = -(1-pwinter)  # to make them look different in colormap

    activeheaders = list(Wheaders) + list(Qheaders)

    # previous year
    df = pd.DataFrame(activeheaders, columns=['Metrics'], index=activeheaders)  # start dataframe
    df['r'] = rwinterp
    df['p'] = pwinterp
    df.drop('Metrics', axis=1, inplace=True)  # get rid of extraneous Metrics column
    # df2 = df.pivot("Metrics", "Previous year", "Previous year")

    # concurrent year
    df2 = pd.DataFrame(activeheaders, columns=['Metrics'], index=activeheaders)  # start dataframe
    df2['r'] = rwinter
    df2['p'] = pwinter
    df2.drop('Metrics', axis=1, inplace=True)  # get rid of extraneous Metrics column

    if whichseason == 'winter':
        cmap = cmocean.cm.waveperiod
    elif whichseason == 'summer':
        cmap = cmocean.cm.waveheight

    # plot previous year calculations
    fig = plt.figure(figsize=(2, 30))
    ax = fig.add_subplot(111)
    sns.heatmap(df, cmap=cmap, ax=ax, vmin=0, vmax=1, cbar=False, xticklabels=['r', 'p'], annot=True)
    # sns.heatmap(np.hstack((rwinterp, -(1-pwinterp))), cmap=cmocean.cm.velocity, ax=ax, vmin=-1, vmax=1, yticklabels=activeheaders, cbar=False)  #, annot=True)
    plt.yticks(rotation=0)
    # being saved on Rainier currently
    fig.savefig('figures/cross/correlations/' + whichseason + 'previous.png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    fig = plt.figure(figsize=(2, 30))
    ax = fig.add_subplot(111)
    sns.heatmap(df2, cmap=cmap, ax=ax, vmin=0, vmax=1, cbar=False, xticklabels=['r', 'p'], annot=True)
    plt.yticks(rotation=0)
    # being saved on Rainier currently
    fig.savefig('figures/cross/correlations/' + whichseason + 'concurrent.png', bbox_inches='tight', dpi=300)
    plt.close(fig)
