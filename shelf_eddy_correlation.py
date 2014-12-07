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
import op
import scipy.stats


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
runcorr = True # run correlation coefficients
doplot = True
explore = False # True to do exploratory plots, False to do polished plot


headers = ('Transport_winter', 'Transport_summer', 
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

if runcorr:
    d = np.loadtxt('table.txt', skiprows=1)

    # check for loop direction problems
    for i in xrange(d.shape[1]):
        if 'Wdm' in headers[i]: # only unwrap if is a mean wind direction metric
            # pdb.set_trace()
            d[:,i] = np.rad2deg(np.unwrap(np.deg2rad(d[:,i])))

    whichseason = 'summer' # 'winter' or 'summer'

    without2008 = False

    if without2008:
        # skip 2008
        d2 = np.vstack((d[0:4,:], d[5:,:]))
        d = d2.copy()

    if whichseason == 'winter':
        T = d[:,0] # the transport is the first column
    elif whichseason == 'summer':
        T = d[:,1] # the transport is the first column
    r = np.empty(d.shape[1]) # to store the correlation coefficients
    p = np.empty(d.shape[1]) # to store the correlation coefficients
    for i in xrange(d.shape[1]):
        comp = d[:,i]
        ind = ~np.isnan(comp)
        r[i], p[i] = scipy.stats.pearsonr(T[ind], comp[ind])

    ind = np.argsort(p) # indices that would sort p in ascending order
    # print the best ones
    for i in xrange(20):
        print headers[ind[i]] + ': p=%1.4f, r=%1.2f' % (p[ind[i]], r[ind[i]])

    # save best entries to plot
    # Top river one
    inan = ~(np.isnan(r)) * ~(r==1)
    ind = np.argmax(r[inan])
    while ('Qbest' not in locals()):
        if 'Q' in np.asarray(headers)[inan][ind]:
            Qbestr = r[inan][ind]
            Qbestp = p[inan][ind]
            Qbest = d[:,inan][:,ind]
            Qbestname = np.asarray(headers)[inan][ind]
        inan *= ~(r==r[inan][ind])
        ind = np.argmax(r[inan])

    # Top wind one
    inan = ~(np.isnan(r)) * ~(r==1)
    ind = np.argmax(r[inan])
    while ('Wbest' not in locals()):
        if 'W' in np.asarray(headers)[inan][ind]:
            Wbestr = r[inan][ind]
            Wbestp = p[inan][ind]
            Wbest = d[:,inan][:,ind]
            Wbestname = np.asarray(headers)[inan][ind]
        inan *= ~(r==r[inan][ind])
        ind = np.argmax(r[inan])


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
                ax.plot(T, Wbest+360, 'o', ms=10, color='0.2', mec='k')
                ax.set_xlim(T.min()-0.1, T.max()+0.1)
                ax.set_xlabel('Transport relative to mean [%]')
                ax.set_ylabel('January-February mean wind direction [deg]')
                ax.text(0.75, 0.15, 'p=%1.4f' % Wbestp, color='r', transform=ax.transAxes, alpha=0.7)
                ax.text(0.75, 0.1, 'r=%1.2f' % Wbestr, color='r', transform=ax.transAxes, alpha=0.7)
                ax.set_frame_on(False)
                fig.savefig('figures/winter-transport-correlations.pdf', bbox_inches='tight')
            elif whichseason == 'summer':
                # with two variables
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111)
                mappable = ax.scatter(d[:,29]/1e5, d[:,117], c=T, cmap='Reds', s=300)
                # ax.set_xlim(T.min()-0.1, T.max()+0.1)
                ax.set_xlabel('March river discharge/10^5 [m$^3\!$/s]')
                ax.set_ylabel('June wind speed variance [m/s]')
                ax.text(0.75, 0.15, '$r^2_{both}$=0.81', color='r', transform=ax.transAxes, alpha=0.7)
                ax.text(0.75, 0.08, '$r^2_{river}$=0.78', color='r', transform=ax.transAxes, alpha=0.7)
                cb = fig.colorbar(mappable)
                cb.set_label('Transport relative to mean [%]')
                ax.set_frame_on(False)
                fig.savefig('figures/summer-transport-correlations2.pdf', bbox_inches='tight')

                # with one variable
                fig = plt.figure(figsize=(6,6))
                ax = fig.add_subplot(111)
                mappable = ax.plot(T, d[:,29]/1e5, 'o', ms=10, color='0.2', mec='k')
                ax.set_xlim(T.min()-0.1, T.max()+0.1)
                ax.set_ylim(3.5, 8.5)
                ax.set_xlabel('Transport relative to mean [%]')
                ax.set_ylabel('March river discharge/10^5 [m$^3\!$/s]')
                ax.text(0.75, 0.15, 'p=%1.4f' % Qbestp, color='r', transform=ax.transAxes, alpha=0.7)
                ax.text(0.75, 0.1, 'r=%1.2f' % Qbestr, color='r', transform=ax.transAxes, alpha=0.7)
                ax.set_frame_on(False)
                fig.savefig('figures/summer-transport-correlations1.pdf', bbox_inches='tight')


if maketable:

    years = np.arange(2004,2015)
    months = np.arange(1,13)

    ## River forcing ##
    # r1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/txla_nesting6/TXLA_river_4dyes_2011.nc')
    r1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012.nc') # use for through 2011
    r2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012_2014.nc') # use for 2012-2014
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
    tr = np.concatenate((tr1, tr2[:]), axis=0)
    datesr = np.concatenate((datesr1, datesr2))
    Q = np.concatenate((Q1, Q2))
    r1.close(); r2.close()

    grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    gridnc = netCDF.Dataset(grid_filename)
    iwind = gridnc.variables['lat_rho'][:]>27.5

    # transport calculated previously
    Tw = np.load('calcs/shelfconn/interannual-wintertransportsum.npz')['transport']
    Ts = np.load('calcs/shelfconn/interannual-summertransportsum.npz')['transport']

    table = np.empty((years.size,17*months.size+2)) # 17 metrics being calculated and also 2 transports

    # Loop through years and write a data point in each column for the row each the year
    for j,year in enumerate(years):
        print year

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

        # Discharge running sum over 1 month, Q1
        Q1 = np.zeros(months.size)
        for i,month in enumerate(months):
            # datetime three months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=1)
            istart = find(datesr<dt)[-1] # starting index
            iend = find(datesr<datetime(year,month,1,0,0,0))[-1] # ending index, current month
            Q1[i] = Q[istart:iend].sum()

        # Discharge running sum over 2 months, Q2
        Q2 = np.zeros(months.size)
        for i,month in enumerate(months):
            # datetime three months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=2)
            istart = find(datesr<dt)[-1] # starting index
            iend = find(datesr<datetime(year,month,1,0,0,0))[-1] # ending index, current month
            Q2[i] = Q[istart:iend].sum()

        # Discharge running sum over three months, Q3
        Q3 = np.zeros(months.size)
        for i,month in enumerate(months):
            # datetime three months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=3)
            istart = find(datesr<dt)[-1] # starting index
            iend = find(datesr<datetime(year,month,1,0,0,0))[-1] # ending index, current month
            Q3[i] = Q[istart:iend].sum()



        ## Wind forcing: averaged over the broad shelf region and taken over 1, 2, and 3 months ##

        # There are multiple file locations
        if year <= 2012:
            w1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year-1) + '.nc')
            w2 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
        elif year == 2013:
            w1 = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_2012.nc')
            w2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
        elif year == 2014:
            w1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
            w2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2014.nc')
            
        # Wind timing
        # pdb.set_trace()
        tw1 = w1.variables['time']
        datesw1 = netCDF.num2date(tw1[:], tw1.units.replace('/', '-'))
        tw2 = w2.variables['time']
        datesw2 = netCDF.num2date(tw2[:], tw2.units.replace('/', '-'))

        # Combine wind info here since files can't be combined in MFDataset
        Uwind1 = w1.variables['Uwind'][:]
        Uwind2 = w2.variables['Uwind'][:]
        Uwind = np.concatenate((Uwind1, Uwind2[1:,:,:]), axis=0)
        Uwind = Uwind[:,iwind] # only use winds from certain area of domain
        del(Uwind1, Uwind2)
        Vwind1 = w1.variables['Vwind'][:]
        Vwind2 = w2.variables['Vwind'][:]
        Vwind = np.concatenate((Vwind1, Vwind2[1:,:,:]), axis=0)
        Vwind = Vwind[:,iwind] # only use winds from certain area of domain
        del(Vwind1, Vwind2)            
        Swind = np.sqrt(Uwind**2 + Vwind**2)

        tw = np.concatenate((tw1, tw2[1:]), axis=0)
        datesw = np.concatenate((datesw1, datesw2[1:]), axis=0)
        w1.close(); w2.close()

        Wsm1 = np.zeros(months.size); Wsv1 = np.zeros(months.size)
        Wdm1 = np.zeros(months.size); Wdv1 = np.zeros(months.size)
        Wsm2 = np.zeros(months.size); Wsv2 = np.zeros(months.size)
        Wdm2 = np.zeros(months.size); Wdv2 = np.zeros(months.size)
        Wsm3 = np.zeros(months.size); Wsv3 = np.zeros(months.size)
        Wdm3 = np.zeros(months.size); Wdv3 = np.zeros(months.size)
        for i, month in enumerate(months):

            # datetime 1 month ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=1)
            istart1 = find(datesw<dt)[-1] # starting index
            # datetime 2 months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=2)
            istart2 = find(datesw<dt)[-1] # starting index
            # datetime 3 months ago
            dt = datetime(year,month,1,0,0,0) - dateutil.relativedelta.relativedelta(months=3)
            istart3 = find(datesw<dt)[-1] # starting index
            iend = find(datesw<datetime(year,month,1,0,0,0))[-1] # ending index, current month

            # Wind magnitude mean, Wsm
            Wsm1[i] = Swind[istart1:iend,:].mean()
            Wsm2[i] = Swind[istart2:iend,:].mean()
            Wsm3[i] = Swind[istart3:iend,:].mean()

            # Wind magnitude var, Wsv
            Wsv1[i] = np.var(Swind[istart1:iend,:])
            Wsv2[i] = np.var(Swind[istart2:iend,:])
            Wsv3[i] = np.var(Swind[istart3:iend,:])

            # Wind direction mean, Wdm
            Wdm1[i] = np.rad2deg(np.arctan2(Vwind[istart1:iend,:].mean(), Uwind[istart1:iend,:].mean()))
            Wdm2[i] = np.rad2deg(np.arctan2(Vwind[istart2:iend,:].mean(), Uwind[istart2:iend,:].mean()))
            Wdm3[i] = np.rad2deg(np.arctan2(Vwind[istart3:iend,:].mean(), Uwind[istart3:iend,:].mean()))

            # Wind direction var, Wdv
            Wdv1[i] = np.var(np.arctan2(Vwind[istart1:iend,:], Uwind[istart1:iend,:]))
            Wdv2[i] = np.var(np.arctan2(Vwind[istart2:iend,:], Uwind[istart2:iend,:]))
            Wdv3[i] = np.var(np.arctan2(Vwind[istart3:iend,:], Uwind[istart3:iend,:]))


        # Combine together into one row (actually put into matrix)
        table[j,:] = np.hstack((Tw[j], Ts[j], Qi, Qcum, Q1, Q2, Q3, Wsm1, Wsm2, Wsm3, 
            Wsv1, Wsv2, Wsv3, Wdm1, Wdm2, Wdm3, Wdv1, Wdv2, Wdv3))

    # write table to file
    with open('table.txt', 'w') as csvfile:
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
# hist(rr, bins=100)