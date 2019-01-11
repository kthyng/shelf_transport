
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
import statsmodels.formula.api as sm
from scipy.signal import argrelextrema

colu = '#218983'  # upcoast color
cold = '#cb6863'  # downcoast color


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
# year = 2004
walong = pd.DataFrame(); wacross = pd.DataFrame()
miss = pd.DataFrame(); atch = pd.DataFrame(); braz = pd.DataFrame()
for year in np.arange(2004,2012):
    # dates = m['ocean_time'].sel(ocean_time=str(year))
    fname_wind = base_wind + str(year) + '.csv'
    fname_river = base_river + str(year) + '.csv'

    walong = pd.concat([walong,pd.read_csv(base_wind + str(year) + '_along.csv', parse_dates=True, index_col=0)], axis=0)
    wacross = pd.concat([wacross,pd.read_csv(base_wind + str(year) + '_across.csv', parse_dates=True, index_col=0)], axis=0)
    miss = pd.concat([miss,pd.read_csv(base_river + str(year) + '_miss.csv', parse_dates=True, index_col=0)], axis=0)
    atch = pd.concat([atch,pd.read_csv(base_river + str(year) + '_atch.csv', parse_dates=True, index_col=0)], axis=0)
    braz = pd.concat([braz,pd.read_csv(base_river + str(year) + '_braz.csv', parse_dates=True, index_col=0)], axis=0)


# since connectivity is by month, average metrics by month too



## Look at spatial correlations
ax = (-sp['2004-01']).T.plot(figsize=(12,4), color='b')
(-ep['2004-01']).T.plot(color='b', linestyle='--',ax=ax)
ax.hlines(0, 0,341, linestyle=':')
(-sn['2004-01']).T.plot(ax=ax, color='orange')
(-en['2004-01']).T.plot(ax=ax, color='orange', linestyle='--')
miss['2004-01'].sum().plot(secondary_y=True, ax=ax)
atch['2004-01'].sum().plot(secondary_y=True, ax=ax)


## look in time for a specific alongcoast location

# find boxes for specific area alongcoast
# miss river plume area
# inds = np.where((dist > 1250) & (dist < 1400))[0]
# sn[inds[10]].plot()
# miss[str(inds[10])].plot(secondary_y=True)

## correlate in time for downcoast connectivity in miss bight and other bumps
# bumpnames = ['miss', 'atch', 'braz', 'main-down', 'main-up']

bumps = {'miss': {'i0': 244, 'i1': 264, 'i2': 290, 'i3': 310},
         'atch': {'i0': 200, 'i1': 220, 'i2': 240, 'i3': 255},
         'braz': {'i0': 131, 'i1': 141, 'i2': 160, 'i3': 175},
         'maindown': {'i0': 70, 'i1': 100, 'i2': 200, 'i3': 220},
         'maindown2': {'i0': 0, 'i1': 30, 'i2': 50, 'i3': 90},
         'mainup': {'i0': 30, 'i1': 60, 'i2': 70, 'i3': 120},
         'mainup2': {'i0': 0, 'i1': 1, 'i2': 30, 'i3': 60}}

mechnames = ['miss','atch','braz','walong','wacross']
mechs = [miss,atch,braz,walong,wacross]
fullnames = []
for bumpname in bumps.keys():
    for mechname in mechnames:
        fullnames.append(bumpname + '_' + str(mechname))  # column names
    fullnames.append(bumpname + '_sum')  # column names

# initialize columns
for fullname in fullnames:
    sn.loc[:,fullname] = np.nan
    sp.loc[:,fullname] = np.nan

# find where break in
for year in np.arange(2004,2012):
    for mon in np.arange(1,13):
        month = '%s-%s' % (year, mon)
        for bump in bumps.keys():
            # bump = 'miss'
            i0, i1, i2, i3 = bumps[bump]['i0'], bumps[bump]['i1'], \
                             bumps[bump]['i2'], bumps[bump]['i3']
            # first  min
            imin0n = sn[month].iloc[:,i0:i1].idxmax(axis=1).values[0]
            # second  min
            imin1n = sn[month].iloc[:,i2:i3].idxmax(axis=1).values[0]
            imin0p = sp[month].iloc[:,i0:i1].idxmin(axis=1).values[0]
            imin1p = sp[month].iloc[:,i2:i3].idxmin(axis=1).values[0]
            # sum over bump
            sn.loc[month,bump + '_sum'] = sn[month].iloc[:,imin0n:imin1n].sum(axis=1)
            sp.loc[month,bump + '_sum'] = sp[month].iloc[:,imin0p:imin1p].sum(axis=1)
            for mech, mechname in zip(mechs,mechnames):
                sn.loc[month,bump + '_' + mechname] = mech[month].iloc[:,imin0n:imin1n].sum(axis=1).sum(axis=0)
                sp.loc[month,bump + '_' + mechname] = mech[month].iloc[:,imin0p:imin1p].sum(axis=1).sum(axis=0)


# # check min locations on the sum by looking at sets of plots
# for year in np.arange(2004,2012):
#     for mon in np.arange(7,9):
#     # mon = 1
#         month = '%s-%s' % (year, mon)
#         fig, ax = plt.subplots(1,1)
#         ax.plot(dist, sp[month].iloc[:,:342].T)
#         for bump in bumps.keys():
#             # bump = 'miss'
#             i0, i1, i2, i3 = bumps[bump]['i0'], bumps[bump]['i1'], \
#                              bumps[bump]['i2'], bumps[bump]['i3']
#             # first  min
#             imin0n = sn[month].iloc[:,i0:i1].idxmax(axis=1).values[0]
#             # second  min
#             imin1n = sn[month].iloc[:,i2:i3].idxmax(axis=1).values[0]
#             imin0p = sp[month].iloc[:,i0:i1].idxmin(axis=1).values[0]
#             imin1p = sp[month].iloc[:,i2:i3].idxmin(axis=1).values[0]
#             print(bump, dist[imin0n], dist[imin1n])
#             # sn[month].iloc[:,imin0n:imin0n+1].plot(marker='o', ax=ax)
#             ax.plot(dist[imin0p], sp[month].iloc[:,imin0p], 'o', alpha=0.7, ms=10)
#             ax.plot(dist[imin1p], sp[month].iloc[:,imin1p], 'x', alpha=0.7, ms=10)
#             ax.set_title(month)

# abs(sn['miss_bump_sum'])['2004':'2007'].plot()
# sn['miss_bump_sum_miss']['2004':'2007'].plot()
# sn['miss_bump_sum_atch']['2004':'2007'].plot()

# # can also then correlate with mississippi river discharge itself
# cols = ['miss_bump_sum','miss_bump_sum_miss','miss_bump_sum_atch',
#         'miss_bump_sum_braz','miss_bump_sum_walong','miss_bump_sum_wacross']
# sn[cols].corr()

for bump in bumps.keys():
    for mechname in mechnames:
        if mechname == bump:
            continue
        name0 = bump + '_sum'
        name1 = bump + '_' + mechname
        result = sm.ols(formula=name0 + " ~ " + name1, data=sn[(sn.index.month==1) | (sn.index.month==2)]).fit()
        if result.pvalues[1] < pcrit:
            print(name0, name1)
            print(result.rsquared)
            # print(result.pvalues[1])
    print('\n')

# # look at correlations
# for col in cols:
#     result = sm.ols(formula="miss_bump_sum ~ " + col, data=sn).fit()
#     if result.pvalues[1] < pcrit:
#         print(col)
#         print(result.rsquared)
#         print(result.pvalues[1])
#
# # best:
# result = sm.ols(formula="miss_bump_sum ~ miss_bump_sum_miss + miss_bump_sum_walong", data=sn).fit()
# print(result.rsquared)
# print(result.pvalues[1:])
#
# MORE HERE
####

## Where along the coastline is the wind correlated in time with connectivity? ##
# include up and downcoast, include along and across wind components
walongmonth = pd.DataFrame(index=sn.index, columns=walong.columns)
wacrossmonth = pd.DataFrame(index=sn.index, columns=wacross.columns)
missmonth = pd.DataFrame(index=sn.index, columns=walong.columns)
atchmonth = pd.DataFrame(index=sn.index, columns=walong.columns)
brazmonth = pd.DataFrame(index=sn.index, columns=walong.columns)
for year in np.arange(2004,2012):
    for mon in np.arange(1,13):
        month = '%s-%s' % (year, mon)
        # sum over month
        walongmonth.loc[month,:] = walong[month].sum(axis=0).values
        wacrossmonth.loc[month,:] = wacross[month].sum(axis=0).values
        missmonth.loc[month,:] = miss[month].sum(axis=0).values
        atchmonth.loc[month,:] = atch[month].sum(axis=0).values
        brazmonth.loc[month,:] = braz[month].sum(axis=0).values

# # df to store correlations along coast
# corrsp = pd.DataFrame(index=walong.columns, columns=['walong','wacross','miss','atch','braz'])
# corrsn = pd.DataFrame(index=walong.columns, columns=['walong','wacross','miss','atch','braz'])
# correp = pd.DataFrame(index=walong.columns, columns=['walong','wacross','miss','atch','braz'])
# corren = pd.DataFrame(index=walong.columns, columns=['walong','wacross','miss','atch','braz'])
# # loop over coastal boxes
# for ibox in range(342):
#     corrsp.loc[str(ibox)] = pd.DataFrame([sp[ibox],walongmonth[str(ibox)],wacrossmonth[str(ibox)],
#                                missmonth[str(ibox)],atchmonth[str(ibox)],brazmonth[str(ibox)]]).T.corr().iloc[0,1:].values
#     corrsn.loc[str(ibox)] = pd.DataFrame([sn[ibox],walongmonth[str(ibox)],wacrossmonth[str(ibox)],
#                                missmonth[str(ibox)],atchmonth[str(ibox)],brazmonth[str(ibox)]]).T.corr().iloc[0,1:].values
#     correp.loc[str(ibox)] = pd.DataFrame([ep[ibox],walongmonth[str(ibox)],wacrossmonth[str(ibox)],
#                                missmonth[str(ibox)],atchmonth[str(ibox)],brazmonth[str(ibox)]]).T.corr().iloc[0,1:].values
#     corren.loc[str(ibox)] = pd.DataFrame([en[ibox],walongmonth[str(ibox)],wacrossmonth[str(ibox)],
#                                missmonth[str(ibox)],atchmonth[str(ibox)],brazmonth[str(ibox)]]).T.corr().iloc[0,1:].values

# abs(corrsp).plot(figsize=(12,4))
# abs(corrsn).plot()

# fig, axes = plt.subplots(2,1,figsize=(12,8))
# axes[0].plot(dist, abs(corrsp['walong']), 'k')
# axes[0].plot(dist, abs(corrsp['wacross']), 'k--')
# axes[0].plot(dist, abs(corrsp['miss']), 'r')
# axes[0].plot(dist, abs(corrsp['atch']), 'r--')
# axes[0].plot(dist, abs(corrsp['braz']), 'r:')
# axes[1].plot(dist, abs(corrsn['walong']), 'k')
# axes[1].plot(dist, abs(corrsn['wacross']), 'k--')
# axes[1].plot(dist, abs(corrsn['miss']), 'r')
# axes[1].plot(dist, abs(corrsn['atch']), 'r--')
# axes[1].plot(dist, abs(corrsn['braz']), 'r:')


# how to separate correlation and causation with river water (which is blown by wind too)
# get associated p values to see where to count
# how to decide when to add over a bump vs. correlate point by point?
cols0 = ['walong','wacross','miss','atch','braz','walong+miss','walong+atch',
         'walong+braz','walong+wacross']
cols1 = ['%s-r2' % col  for col in cols0]
cols2 = ['%s-p' % col  for col in cols0]
cols3 = ['%s-p2' % col  for col in cols0 if '+' in col]
cols = cols1 + cols2 + cols3

corrspall = pd.DataFrame(index=walong.columns, columns=cols)
corrsnall = pd.DataFrame(index=walong.columns, columns=cols)
# loop over coastal boxes
for ibox in range(342):
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walongmonth[str(ibox)],
                             'wacross': wacrossmonth[str(ibox)],
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walongmonth[str(ibox)],
                             'wacross': wacrossmonth[str(ibox)],
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    for col in cols0:
        result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
        corrspall.loc[str(ibox),col+'-r2'] = result.rsquared_adj
        corrspall.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrspall.loc[str(ibox),col+'-p2'] = result.pvalues[2]
        result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
        corrsnall.loc[str(ibox),col+'-r2'] = result.rsquared_adj
        corrsnall.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrsnall.loc[str(ibox),col+'-p2'] = result.pvalues[2]
# print(result.summary())

# only show r^2 if p is low enough, p<0.01
pcrit = 0.01
fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrspall[col+'-p']<pcrit) & (corrspall[col+'-p2']<pcrit)
    else:
        inds = corrspall[col+'-p']<pcrit
    x = dist.copy()
    y = corrspall[col+'-r2'].values
    y[~inds] = np.nan
    ax.plot(x, y, lw=2, label=col)
fig.legend()
ax.set_xlabel('Along-coast distance [km]')
ax.set_ylabel('r$^2$')
# also plot overall connectivity
ax2 = ax.twinx()
ax2.fill_between(dist, sp.iloc[:,:342].sum(axis=0), color=colu, alpha=0.3)
fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_upcoast.png', bbox_inches='tight')
fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_upcoast.pdf', bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrsnall[col+'-p']<pcrit) & (corrsnall[col+'-p2']<pcrit)
    else:
        inds = corrsnall[col+'-p']<pcrit
    x = dist.copy()
    y = corrsnall[col+'-r2'].values
    y[~inds] = np.nan
    ax.plot(x, y, lw=2, label=col)
fig.legend()
ax.set_xlabel('Along-coast distance [km]')
ax.set_ylabel('r$^2$')
# also plot overall connectivity
ax2 = ax.twinx()
ax2.fill_between(dist, abs(sn.iloc[:,:342].sum(axis=0)), color=cold, alpha=0.3)
fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_downcoast.png', bbox_inches='tight')
fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_downcoast.pdf', bbox_inches='tight')

# corrspall[['walong-r2','miss-r2','walong+miss-r2']].plot()

# seem to be able to pull out brazos, atch, convergence in cross shelf wind,
# but not miss — seem to need to add over bump for that
# how to decide when need to add over bump?
# maybe can add over all bumps and look for good statistical values there too

# do with more metrics in time included from rainier

# fig, axes = plt.subplots(2,1,figsize=(12,8))
# corrspall[['walong-r2','miss-r2','walong+miss-r2']].plot(ax=axes[0])
# corrspall[['walong-p','miss-p','walong+miss-p','walong+miss-p2']].plot(ax=axes[1])

# also look at interannual variability

# separate out seasons
iwin = (sp.index.month == 1) | (sp.index.month == 2)
isum = (sp.index.month == 7) | (sp.index.month == 8)
corrspwin = pd.DataFrame(index=walong.columns, columns=cols)
corrsnwin = pd.DataFrame(index=walong.columns, columns=cols)
corrspsum = pd.DataFrame(index=walong.columns, columns=cols)
corrsnsum = pd.DataFrame(index=walong.columns, columns=cols)
# loop over coastal boxes
for ibox in range(342):
    dftpwin = pd.DataFrame(data={'sp': sp[ibox][iwin],'walong': walongmonth[str(ibox)][iwin],
                             'wacross': wacrossmonth[str(ibox)][iwin],
                             'miss': missmonth[str(ibox)][iwin],
                             'atch':atchmonth[str(ibox)][iwin],
                             'braz':brazmonth[str(ibox)][iwin]}, dtype=float64)
    dftnwin = pd.DataFrame(data={'sn': sn[ibox][iwin],'walong': walongmonth[str(ibox)][iwin],
                             'wacross': wacrossmonth[str(ibox)][iwin],
                             'miss': missmonth[str(ibox)][iwin],
                             'atch':atchmonth[str(ibox)][iwin],
                             'braz':brazmonth[str(ibox)][iwin]}, dtype=float64)
    dftpsum = pd.DataFrame(data={'sp': sp[ibox][isum],'walong': walongmonth[str(ibox)][isum],
                             'wacross': wacrossmonth[str(ibox)][isum],
                             'miss': missmonth[str(ibox)][isum],
                             'atch':atchmonth[str(ibox)][isum],
                             'braz':brazmonth[str(ibox)][isum]}, dtype=float64)
    dftnsum = pd.DataFrame(data={'sn': sn[ibox][isum],'walong': walongmonth[str(ibox)][isum],
                             'wacross': wacrossmonth[str(ibox)][isum],
                             'miss': missmonth[str(ibox)][isum],
                             'atch':atchmonth[str(ibox)][isum],
                             'braz':brazmonth[str(ibox)][isum]}, dtype=float64)
    for col in cols0:
        result = sm.ols(formula="sp ~ " + col, data=dftpwin).fit()
        corrspwin.loc[str(ibox),col+'-r2'] = result.rsquared
        corrspwin.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrspwin.loc[str(ibox),col+'-p2'] = result.pvalues[2]

        result = sm.ols(formula="sn ~ " + col, data=dftnwin).fit()
        corrsnwin.loc[str(ibox),col+'-r2'] = result.rsquared
        corrsnwin.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrsnwin.loc[str(ibox),col+'-p2'] = result.pvalues[2]

        result = sm.ols(formula="sp ~ " + col, data=dftpsum).fit()
        corrspsum.loc[str(ibox),col+'-r2'] = result.rsquared
        corrspsum.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrspsum.loc[str(ibox),col+'-p2'] = result.pvalues[2]

        result = sm.ols(formula="sn ~ " + col, data=dftnsum).fit()
        corrsnsum.loc[str(ibox),col+'-r2'] = result.rsquared
        corrsnsum.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrsnsum.loc[str(ibox),col+'-p2'] = result.pvalues[2]
# print(result.summary())

# only show r^2 if p is low enough, p<0.01
pcrit = 0.01
fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrspwin[col+'-p']<pcrit) & (corrspwin[col+'-p2']<pcrit)
        ax.plot(dist[inds], corrspwin[col+'-r2'][inds], marker='.')
    else:
        inds = corrspwin[col+'-p']<pcrit
        ax.plot(dist[inds], corrspwin[col+'-r2'][inds], marker='.')
fig.legend()

fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrsnwin[col+'-p']<pcrit) & (corrsnwin[col+'-p2']<pcrit)
        ax.plot(dist[inds], corrsnwin[col+'-r2'][inds], marker='.')
    else:
        inds = corrsnwin[col+'-p']<pcrit
        ax.plot(dist[inds], corrsnwin[col+'-r2'][inds], marker='.')
fig.legend()

fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrspsum[col+'-p']<pcrit) & (corrspsum[col+'-p2']<pcrit)
        ax.plot(dist[inds], corrspsum[col+'-r2'][inds], marker='.')
    else:
        inds = corrspsum[col+'-p']<pcrit
        ax.plot(dist[inds], corrspsum[col+'-r2'][inds], marker='.')
fig.legend()

fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrsnsum[col+'-p']<pcrit) & (corrsnsum[col+'-p2']<pcrit)
        ax.plot(dist[inds], corrsnsum[col+'-r2'][inds], marker='.')
    else:
        inds = corrsnsum[col+'-p']<pcrit
        ax.plot(dist[inds], corrsnsum[col+'-r2'][inds], marker='.')
fig.legend()



## correlate based on between bumps ##

# seems to hard to generically find the bumps from the local maxima/minima -
# too much variability
# I think I need to specify windows to look in for known bumps

# downcoast, winter
# TEST ONE BOX BY SEASON AND DIRECTION
iboxes = np.where((dist > 300) & (dist < 500))[0]  # which boxes
y = []
x0 = []; x1 = []
for year in np.arange(2004,2012):
    for month in np.arange(1,13):
        mon = '%s-%s' % (year, month)
        ibump0 = argrelextrema(sn[mon].iloc[:,iboxes].values, np.greater_equal, axis=1, order=20)[1]
        plt.plot(dist, sn[mon].iloc[:,:342].T)
        plt.plot(dist[iboxes[ibump0]], sn[mon].iloc[:,iboxes[ibump0]].T, 'o')
        y.append(sn[mon].iloc[:,:iboxes[ibump0][0]].sum(axis=1))
        x0.append(walong[mon].iloc[:,:iboxes[ibump0][0]].sum(axis=1).sum(axis=0))
        x1.append(wacross[mon].iloc[:,:iboxes[ibump0][0]].sum(axis=1).sum(axis=0))

result = sm.ols(formula="y ~ x0", data={'y': y, 'x0': x0, 'x1': x1}).fit()
print(result.rsquared)
print(result.pvalues[1:])
result = sm.ols(formula="y ~ x1", data={'y': y, 'x0': x0, 'x1': x1}).fit()
print(result.rsquared)
print(result.pvalues[1:])
result = sm.ols(formula="y ~ x0+x1", data={'y': y, 'x0': x0, 'x1': x1}).fit()
print(result.rsquared)
print(result.pvalues[1:])




# find local max/minima in connectivity
n = 20 # number of points to be checked before and after
imin = argrelextrema(sp['2004-07'].iloc[:,:342].values, np.less_equal, axis=1, order=n)[1]
plt.plot(dist, sp['2004-07'].iloc[:,:342].T)
plt.plot(dist[imin], sp['2004-07'].iloc[:,imin].T, 'o')

n = 20
figure()
imax = argrelextrema(sn['2004-01'].iloc[:,:342].values, np.greater_equal, axis=1, order=n)[1]
plt.plot(dist, sn['2004-01'].iloc[:,:342].T)
plt.plot(dist[imax], sn['2004-01'].iloc[:,imax].T, 'o')

imax = argrelextrema(sn['2005-01'].iloc[:,:342].values, np.greater_equal, axis=1, order=n)[1]
plt.plot(dist, sn['2005-01'].iloc[:,:342].T)
plt.plot(dist[imax], sn['2005-01'].iloc[:,imax].T, 'o')

imax = argrelextrema(sn['2006-01'].iloc[:,:342].values, np.greater_equal, axis=1, order=n)[1]
plt.plot(dist, sn['2006-01'].iloc[:,:342].T)
plt.plot(dist[imax], sn['2006-01'].iloc[:,imax].T, 'o')
####



## try correlating by a set number of group boxes ##
cols0 = ['walong','wacross','miss','atch','braz','walong+miss','walong+atch',
         'walong+braz','walong+wacross']
cols1 = ['%s-r2' % col  for col in cols0]
cols2 = ['%s-p' % col  for col in cols0]
cols3 = ['%s-p2' % col  for col in cols0 if '+' in col]
cols = cols1 + cols2 + cols3
dd = 20

corrspall = pd.DataFrame(index=walong.columns[::dd], columns=cols)
corrsnall = pd.DataFrame(index=walong.columns[::dd], columns=cols)
# loop over coastal boxes
for ibox in range(0,342-dd,dd):
    dftp = pd.DataFrame(data={'sp': sp.loc[:,ibox:ibox+dd-1].sum(axis=1, skipna=False),'walong': walongmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'wacross': wacrossmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'miss': missmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'atch':atchmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'braz':brazmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False)}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn.iloc[:,ibox:ibox+dd-1].sum(axis=1, skipna=False),'walong': walongmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'wacross': wacrossmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'miss': missmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'atch':atchmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False),
                             'braz':brazmonth.loc[:,str(ibox):str(ibox+dd-1)].sum(axis=1, skipna=False)}, dtype=float64)
    for col in cols0:
        result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
        corrspall.loc[str(ibox),col+'-r2'] = result.rsquared
        corrspall.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrspall.loc[str(ibox),col+'-p2'] = result.pvalues[2]
        result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
        corrsnall.loc[str(ibox),col+'-r2'] = result.rsquared
        corrsnall.loc[str(ibox),col+'-p'] = result.pvalues[1]
        if '+' in col:
            corrsnall.loc[str(ibox),col+'-p2'] = result.pvalues[2]
# print(result.summary())



# only show r^2 if p is low enough, p<0.01
pcrit = 0.01
fig, ax = plt.subplots(1,1,figsize=(12,4))
for col in cols0:
    if '+' in col:
        inds = (corrspall[col+'-p']<pcrit) & (corrspall[col+'-p2']<pcrit)
        # if inds.diff().sum() > 5:
        ax.plot(dist[::dd][inds], corrspall[col+'-r2'][inds], marker='.')
    else:
        inds = corrspall[col+'-p']<pcrit
        ax.plot(dist[::dd][inds], corrspall[col+'-r2'][inds], marker='.')
fig.legend()

fig, ax = plt.subplots(1,1,figsize=(12,4))
# col = cols0[0]
for col in cols0:
    if '+' in col:
        inds = (corrsnall[col+'-p']<pcrit) & (corrsnall[col+'-p2']<pcrit)
        # if inds.diff().sum() > 5:
        ax.plot(dist[::dd][inds], corrsnall[col+'-r2'][inds], marker='.')
    else:
        inds = corrsnall[col+'-p']<pcrit
        ax.plot(dist[::dd][inds], corrsnall[col+'-r2'][inds], marker='.')
fig.legend()
####
