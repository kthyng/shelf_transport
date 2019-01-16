
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
        # mean over month
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


# # how often is the wind along/across coast and in what directions
# fig = plt.figure(figsize=(7, 4))
# ax = fig.add_axes(params[loc0]['axes'], projection=merc)
# ax.set_frame_on(False) # kind of like it without the box
# ax.set_extent([-98, -88.5, 25.5, 30.3], pc)
# gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
# # the following two make the labels look like lat/lon format
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.xlabels_bottom = False  # turn off labels where you don't want them
# gl.ylabels_right = False
# gl.xlabels_top = False  # turn off labels where you don't want them
# gl.ylabels_left = False
# ax.add_feature(land_10m, facecolor='0.85')
# # ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
# ax.add_feature(cartopy.feature.BORDERS, facecolor='0.8')
# ax.add_feature(states_provinces, edgecolor='gray')





# how to separate correlation and causation with river water (which is blown by wind too)
# get associated p values to see where to count
# how to decide when to add over a bump vs. correlate point by point?
cols0 = ['walong','wacross','miss','atch','braz','walong+miss','walong+atch',
         'walong+braz','walong+wacross']
cols1 = ['%s-r2' % col  for col in cols0]
cols2 = ['%s-p' % col  for col in cols0]
cols3 = ['%s-p2' % col  for col in cols0 if '+' in col]
cols = cols1 + cols2 + cols3

# loop over range of wind angle shifts to see impact on correlations alongcoast
for rot in np.arange(0,360,15):  # every 15 deg

    corrspall = pd.DataFrame(index=walong.columns, columns=cols)
    corrsnall = pd.DataFrame(index=walong.columns, columns=cols)
    # loop over coastal boxes
    for ibox in range(342):
        # rotate wind components by angle rot
        walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
                    - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
        wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
                    + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
        dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                                 'wacross': wacrossuse,
                                 'miss': missmonth[str(ibox)],
                                 'atch':atchmonth[str(ibox)],
                                 'braz':brazmonth[str(ibox)]}, dtype=float64)
        dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                                 'wacross': wacrossuse,
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
    fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_upcoast_rot%i.png' % rot, bbox_inches='tight')
    fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_upcoast_rot%i.pdf' % rot, bbox_inches='tight')

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
    fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_downcoast_rot%i.png' % rot, bbox_inches='tight')
    fig.savefig('figures/alongcoastconn/alongcoast_r2_adj_downcoast_rot%i.pdf' % rot, bbox_inches='tight')

# corrspall[['walong-r2','miss-r2','walong+miss-r2']].plot()


# look at scatter plot of a particular coastal box to see what is being correlated
# ibox = 75
base = 'figures/alongcoastconn/test/windvec'
os.makedirs(base, exist_ok=True)
# for rot in np.arange(0,360,15):  # every 15 deg
rot = 0
for ibox in [20,75,150]:
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    vmax = max((abs(dftn['sn']).max(), abs(dftp['sp']).max()))

    fig, axes = plt.subplots(1,2)#, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.09, bottom=0.22, right=0.97, top=0.94, wspace=0.11)
    axes[0].set_ylabel('wind: across')
    axes[0].set_xlabel('wind: along')
    axes[1].set_xlabel('wind: along')
    axes[0].scatter(dftn['walong'], dftn['wacross'], c=dftn['sn'], s=100,cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
    mappable = axes[1].scatter(dftp['walong'], dftp['wacross'], c=dftp['sp'], s=100,cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
    axes[0].hlines(0, dftp['walong'].min(), dftp['walong'].max(), linestyle=':')
    axes[0].vlines(0, dftp['wacross'].min(), dftp['wacross'].max(), linestyle=':')
    axes[1].hlines(0, dftp['walong'].min(), dftp['walong'].max(), linestyle=':')
    axes[1].vlines(0, dftp['wacross'].min(), dftp['wacross'].max(), linestyle=':')
    axes[0].set_ylim(-3,3); axes[0].set_xlim(-3,3);
    axes[1].set_ylim(-3,3); axes[1].set_xlim(-3,3);
    # axes[0].axis('equal'); axes[1].axis('equal')
    cax = fig.add_axes([0.2, 0.09, 0.7, 0.025])
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Connectivity')
    # add wind stats to plot
    for i, col in enumerate(['walong', 'wacross', 'walong+wacross']):
        result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
        ax = axes[1]
        ax.text(0.01, 0.95-i*0.2, '%s' % col, transform=ax.transAxes)
        ax.text(0.01, 0.9-i*0.2, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
        # print('r$^2$=%0.2f' % result.rsquared_adj)
        ax.text(0.01, 0.85-i*0.2, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
        if '+' in col:
            ax.text(0.01, 0.8-i*0.2, 'p2=%0.2f' % result.pvalues[2], transform=ax.transAxes)

        result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
        ax = axes[0]
        ax.text(0.01, 0.95-i*0.2, '%s' % col, transform=ax.transAxes)
        ax.text(0.01, 0.9-i*0.2, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
        # print('r$^2$=%0.2f' % result.rsquared_adj)
        ax.text(0.01, 0.85-i*0.2, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
        if '+' in col:
            ax.text(0.01, 0.8-i*0.2, 'p2=%0.2f' % result.pvalues[2], transform=ax.transAxes)
    fig.suptitle('ibox=%2d, dist=%4d' % (ibox,dist[ibox]))
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)


# separate by up/downcoast dominant first for two time series
ibox = 75
base = 'figures/alongcoastconn/test/windvec/sep'
os.makedirs(base, exist_ok=True)
for rot in np.arange(0,360,15):  # every 15 deg
# rot = 0
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    vmax = max((abs(dftn['sn']).max(), abs(dftp['sp']).max()))

    iup = abs(dftp['sp']) > abs(dftn['sn'])
    idown = ~iup  # time indices when downcoast dominant

    fig, axes = plt.subplots(1,2, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.09, bottom=0.22, right=0.97, top=0.94, wspace=0.11)
    axes[0].set_ylabel('wind: across')
    axes[0].set_xlabel('wind: along')
    axes[1].set_xlabel('wind: along')
    axes[0].scatter(dftn['walong'][idown], dftn['wacross'][idown], c=dftn['sn'][idown], s=100,cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
    mappable = axes[1].scatter(dftp['walong'][iup], dftp['wacross'][iup], c=dftp['sp'][iup], s=100,cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
    axes[0].hlines(0, dftp['walong'][idown].min(), dftp['walong'][idown].max(), linestyle=':')
    axes[0].vlines(0, dftp['wacross'][iup].min(), dftp['wacross'][iup].max(), linestyle=':')
    axes[1].hlines(0, dftp['walong'][idown].min(), dftp['walong'][idown].max(), linestyle=':')
    axes[1].vlines(0, dftp['wacross'][iup].min(), dftp['wacross'][iup].max(), linestyle=':')
    cax = fig.add_axes([0.2, 0.09, 0.7, 0.025])
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')
    cb.set_label('Connectivity')
    # add wind stats to plot
    for i, col in enumerate(['walong', 'wacross', 'walong+wacross']):
        result = sm.ols(formula="sp ~ " + col, data=dftp[iup]).fit()
        ax = axes[1]
        ax.text(0.01, 0.95-i*0.2, '%s' % col, transform=ax.transAxes)
        ax.text(0.01, 0.9-i*0.2, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
        # print('r$^2$=%0.2f' % result.rsquared_adj)
        ax.text(0.01, 0.85-i*0.2, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
        if '+' in col:
            ax.text(0.01, 0.8-i*0.2, 'p2=%0.2f' % result.pvalues[2], transform=ax.transAxes)

        result = sm.ols(formula="sn ~ " + col, data=dftn[idown]).fit()
        ax = axes[0]
        ax.text(0.01, 0.95-i*0.2, '%s' % col, transform=ax.transAxes)
        ax.text(0.01, 0.9-i*0.2, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
        # print('r$^2$=%0.2f' % result.rsquared_adj)
        ax.text(0.01, 0.85-i*0.2, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
        if '+' in col:
            ax.text(0.01, 0.8-i*0.2, 'p2=%0.2f' % result.pvalues[2], transform=ax.transAxes)
    fig.suptitle('rotation angle=%2d' % rot)
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)






# several examples follow, which I am now using to systematically find the angle
# for each coastal box for which 1. wacross>0 and 2. r^2 is minimized.
# these tend to give the best relationship for walong once all are

# another scatter plot attempt
ibox = 120
base = 'figures/alongcoastconn/test/'
os.makedirs(base, exist_ok=True)
# for rot in np.arange(0,360,15):  # every 15 deg
rot = 0
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)


    fig, axes = plt.subplots(2,2, sharex=True, sharey=True)
    vmax = dftp['atch'].max(); vmin = dftp['atch'].min()
    for i, row in enumerate([dftn['sn'], dftp['sp']]):
        for j, col in enumerate(['walong', 'wacross']):
            ax = axes[i,j]
            ax.scatter(dftp[col], row, c=dftp['atch'], s=50, cmap='cmo.amp', vmin=vmin, vmax=vmax)#, norm=matplotlib.colors.LogNorm())
            ax.hlines(0, dftp[col].min(), dftp[col].max())
            ax.vlines(0, row.min(), row.max())
            if i==0:
                result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
            elif i==1:
                result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)

    axes[0,0].set_ylabel('downcoast')
    axes[1,0].set_ylabel('upcoast')
    axes[1,0].set_xlabel('walong')
    axes[1,1].set_xlabel('wacross')
    fig.suptitle('rotation angle=%2d' % rot)
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)

# another scatter plot attempt
ibox = 150
base = 'figures/alongcoastconn/test/'
os.makedirs(base, exist_ok=True)
for rot in np.arange(0,360,15):  # every 15 deg
# rot = 45
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)


    fig, axes = plt.subplots(2,2, sharex=True, sharey=True)
    vmax = dftp['atch'].max(); vmin = dftp['atch'].min()
    for i, row in enumerate([dftn['sn'], dftp['sp']]):
        for j, col in enumerate(['walong', 'wacross']):
            ax = axes[i,j]
            ax.scatter(dftp[col], row, c=dftp['atch'], s=50, cmap='cmo.amp', vmin=vmin, vmax=vmax)#, norm=matplotlib.colors.LogNorm())
            ax.hlines(0, dftp[col].min(), dftp[col].max())
            ax.vlines(0, row.min(), row.max())
            if i==0:
                result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
            elif i==1:
                result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)

    axes[0,0].set_ylabel('downcoast')
    axes[1,0].set_ylabel('upcoast')
    axes[1,0].set_xlabel('walong')
    axes[1,1].set_xlabel('wacross')
    fig.suptitle('rotation angle=%2d' % rot)
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)

# 315 angle shows only connectivity when onshore, which makes sense
# check a box between 250 and 500 to see how that looks
ibox = 75
base = 'figures/alongcoastconn/test/'
os.makedirs(base, exist_ok=True)
for rot in np.arange(0,360,15):  # every 15 deg
# rot = 45
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)


    fig, axes = plt.subplots(2,2, sharex=True, sharey=True)
    # i = 0
    for i, row in enumerate([dftn['sn'], dftp['sp']]):
        for j, col in enumerate(['walong', 'wacross']):
            ax = axes[i,j]
            ax.plot(dftp[col], row, 'o')
            ax.hlines(0, dftp[col].min(), dftp[col].max())
            ax.vlines(0, row.min(), row.max())
            if i==0:
                result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
            elif i==1:
                result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)

    axes[0,0].set_ylabel('downcoast')
    axes[1,0].set_ylabel('upcoast')
    axes[1,0].set_xlabel('walong')
    axes[1,1].set_xlabel('wacross')
    fig.suptitle('rotation angle=%2d' % rot)
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)


ibox = 20
base = 'figures/alongcoastconn/test/'
os.makedirs(base, exist_ok=True)
for rot in np.arange(0,360,15):  # every 15 deg
# rot = 45
    walonguse = walongmonth[str(ibox)]*np.cos(np.deg2rad(rot)) \
            - wacrossmonth[str(ibox)]*np.sin(np.deg2rad(rot))
    wacrossuse = walongmonth[str(ibox)]*np.sin(np.deg2rad(rot)) \
            + wacrossmonth[str(ibox)]*np.cos(np.deg2rad(rot))
    dftp = pd.DataFrame(data={'sp': sp[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)
    dftn = pd.DataFrame(data={'sn': sn[ibox],'walong': walonguse,
                             'wacross': wacrossuse,
                             'miss': missmonth[str(ibox)],
                             'atch':atchmonth[str(ibox)],
                             'braz':brazmonth[str(ibox)]}, dtype=float64)


    fig, axes = plt.subplots(2,2, sharex=True, sharey=True)
    # i = 0
    for i, row in enumerate([dftn['sn'], dftp['sp']]):
        for j, col in enumerate(['walong', 'wacross']):
            ax = axes[i,j]
            ax.plot(dftp[col], row, 'o')
            ax.hlines(0, dftp[col].min(), dftp[col].max())
            ax.vlines(0, row.min(), row.max())
            if i==0:
                result = sm.ols(formula="sn ~ " + col, data=dftn).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)
            elif i==1:
                result = sm.ols(formula="sp ~ " + col, data=dftp).fit()
                ax.text(0.01, 0.9, '%s' % col, transform=ax.transAxes)
                ax.text(0.01, 0.8, 'r$^2$=%0.2f' % result.rsquared_adj, transform=ax.transAxes)
                ax.text(0.01, 0.7, 'p=%0.2f' % result.pvalues[1], transform=ax.transAxes)

    axes[0,0].set_ylabel('downcoast')
    axes[1,0].set_ylabel('upcoast')
    axes[1,0].set_xlabel('walong')
    axes[1,1].set_xlabel('wacross')
    fig.suptitle('rotation angle=%2d' % rot)
    fig.savefig('%s/ibox_%i_rot_%2d' % (base, ibox, rot))
    plt.close(fig)

# should i still keep up and downcoast separately or how to deal otherwise?
# Maybe I need to account for time lag of wind, or how long it has been blowing



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




# correlate with high time res
# this didn't really work out
walong2 = walong.resample('4H').interpolate()
wacross2 = wacross.resample('4H').interpolate()

d = np.load('calcs/alongcoastconn/conn_in_time/lines_by_start.npz')
sp2 = d['startingp']; sn2 = d['startingn']
# there are a bunch of -0.0625 in downcoast sn. Set to nan because I think they
# are a replacement for 0 and for high time res they don't look right to include
sn2[sn2 == -0.0625] = np.nan
# run correlations for simulations across season
r2 = []
vmax = abs(np.nanmin(sn2[:11,:,ibox]))
for dstart in pd.date_range(start=dstartoverall, end=pd.Timestamp('2004-1-11'), freq= '1D'):
    # dstart = pd.Timestamp('2004-01-02'); dt = pd.Timedelta('30 days')
    dstartoverall = pd.Timestamp('2004-01-01')
    dend = dstart + dt - pd.Timedelta('4 hours')
    ibox = 150
    idays = (dstart - dstartoverall).days  # days after 2st possible date, which gives index

    fig = plt.figure()
    x = walong2.loc[dstart:dend,str(ibox)]
    y = wacross2.loc[dstart:dend,str(ibox)]
    mappable = plt.scatter(x, y, c=sn2[idays,:,ibox], s=100, cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
    # plt.scatter(walong2.loc[dstart:dend,str(ibox)], wacross2.loc[dstart:dend,str(ibox)], c=sp2[idays,:,ibox], s=100, cmap='cmo.curl_r')
    fig.colorbar(mappable)
    plt.hlines(0, x.min(), x.max(), linestyle=':')
    plt.vlines(0, y.min(), y.max(), linestyle=':')
    plt.ylabel('wacross')
    plt.xlabel('walong')

    result = sm.ols(formula="sn ~ walong", data={'walong': walong2.loc[dstart:dend,str(ibox)], 'sn': sn2[idays,:,ibox]}).fit()
    r2.append(result.rsquared)
    # print(result.rsquared)
    # print(result.pvalues[1:])



# try with time res of 1 per simulation
walong2 = walong.resample('1D').interpolate()
wacross2 = wacross.resample('1D').interpolate()
walong3 = walong.resample('1MS').sum()
wacross3 = wacross.resample('1MS').sum()

d = np.load('calcs/alongcoastconn/conn_in_time/lines_by_start_sim.npz')
sp2 = d['startingp']; sn2 = d['startingn']
index = pd.DatetimeIndex(d['dates'])
sn3 = pd.DataFrame(index=index, columns=np.arange(0,len(dist)), data=sn2)
sp3 = pd.DataFrame(index=index, columns=np.arange(0,len(dist)), data=sp2)

sp4 = sp3.resample('1MS').mean()
sn4 = sn3.resample('1MS').mean()

dend = pd.Timestamp('2011-12-31')
ibox = 75
# idays = (dstart - dstartoverall).days  # days after 2st possible date, which gives index
iday = np.where(index == dend)[0][0]

which = 'try3'  # 'orig', 'try1'
fig = plt.figure()
# LOOKS LIKE I CAN GET EVEN BETTER RESULTS NOW, BUT WHY ISN'T DAILY BETTER THAN
# MONTHLY?
# original
if which == 'orig':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn[ibox]
    title = 'original'
# all new, by day
elif which == 'try1':
    x = walong2.loc[:,str(ibox)]
    y = wacross2.loc[:,str(ibox)]
    z = sn3[:walong2.index[-1]][ibox]
    title = 'all new, by day'
# new daily, resampling to months
elif which == 'try2':
    x = walong3.loc[:,str(ibox)]
    y = wacross3.loc[:,str(ibox)]
    z = sn4[:walong3.index[-1]][ibox]
    title = 'new daily, resampling to monthly'
# new daily, using previous wind monthly sum
elif which == 'try3':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn4[ibox]
    title = 'new daily, previous wind sum'
#
vmax = abs(np.nanmin(z))

mappable = plt.scatter(x, y, c=z, s=100, cmap='cmo.curl_r', vmin=-vmax, vmax=vmax);
# mappable = plt.scatter(x, y, c=sn2[:iday+1,ibox], s=100, cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
# plt.scatter(walong2.loc[dstart:dend,str(ibox)], wacross2.loc[dstart:dend,str(ibox)], c=sp2[idays,:,ibox], s=100, cmap='cmo.curl_r')
fig.colorbar(mappable)
plt.hlines(0, x.min(), x.max(), linestyle=':')
plt.vlines(0, y.min(), y.max(), linestyle=':')
plt.ylabel('wacross')
plt.xlabel('walong')

# result = sm.ols(formula="sn ~ walong", data=dftn).fit()
dftemp = pd.DataFrame(index=x.index, data={'sn': z, 'walong': x, 'wacross': y}, dtype=float)
result = sm.ols(formula="sn ~ walong", data=dftemp).fit()
# r2.append(result.rsquared)
print(result.rsquared)
# print(result.pvalues[1:])
plt.title(title + ', r$^2$=%0.2f' % result.rsquared)




d = np.load('calcs/alongcoastconn/conn_in_time/lines_by_start_4hours.npz')
# number of simulations (matching dates) x time in simulation (4 hours) x alongcoast conn box
sp5 = d['startingp']; sn5 = d['startingn']
# change to monthly to match sp/sn:
# leave 0 axis, sum of axis 1, leave axis 2
# sum across month
sp6 = pd.DataFrame(index=d['dates'], data=sp5.sum(axis=1))
sp6 = sp6.resample('1MS').sum()
sn6 = pd.DataFrame(index=d['dates'], data=sn5.sum(axis=1))
sn6 = sn6.resample('1MS').sum()
# change to monthly to match sp/sn:
# leave 0 axis, sum of axis 1, leave axis 2
# mean across month
sp7 = pd.DataFrame(index=d['dates'], data=sp5.sum(axis=1))
sp7 = sp7.resample('1MS').mean()
sn7 = pd.DataFrame(index=d['dates'], data=sn5.sum(axis=1))
sn7 = sn7.resample('1MS').mean()
# change to monthly to match sp/sn:
# leave 0 axis, mean of axis 1, leave axis 2
# mean across month
sp8 = pd.DataFrame(index=d['dates'], data=sp5.mean(axis=1))
sp8 = sp8.resample('1MS').mean()
sn8 = pd.DataFrame(index=d['dates'], data=sn5.mean(axis=1))
sn8 = sn8.resample('1MS').mean()
# change to monthly to match sp/sn:
# leave 0 axis, mean of axis 1, leave axis 2
# mean across month
sp9 = pd.DataFrame(index=d['dates'], data=sp5.mean(axis=1))
sp9 = sp9.resample('1MS').sum()
sn9 = pd.DataFrame(index=d['dates'], data=sn5.mean(axis=1))
sn9 = sn9.resample('1MS').sum()

which = 'try6'
# try to match original plot to make sure doing same process
if which == 'orig':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn[ibox]
    title = 'original'
elif which == 'try4':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn6[ibox]
    title = 'trying to match, with sum/sum'
elif which == 'try5':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn7[ibox]
    title = 'trying to match, with sum/mean'
elif which == 'try6':
    x = walongmonth.loc[:,str(ibox)]
    y = wacrossmonth.loc[:,str(ibox)]
    z = sn8[ibox]
    title = 'trying to match, with mean/mean'

vmax = abs(np.nanmin(z))

fig = plt.figure()
mappable = plt.scatter(x, y, c=z, s=100, cmap='cmo.curl_r', vmin=-vmax, vmax=vmax);
# mappable = plt.scatter(x, y, c=sn2[:iday+1,ibox], s=100, cmap='cmo.curl_r', vmin=-vmax, vmax=vmax)
# plt.scatter(walong2.loc[dstart:dend,str(ibox)], wacross2.loc[dstart:dend,str(ibox)], c=sp2[idays,:,ibox], s=100, cmap='cmo.curl_r')
fig.colorbar(mappable)
plt.hlines(0, x.min(), x.max(), linestyle=':')
plt.vlines(0, y.min(), y.max(), linestyle=':')
plt.ylabel('wacross')
plt.xlabel('walong')

# result = sm.ols(formula="sn ~ walong", data=dftn).fit()
dftemp = pd.DataFrame(index=x.index, data={'sn': z, 'walong': x, 'wacross': y}, dtype=float)
result = sm.ols(formula="sn ~ walong", data=dftemp).fit()
# r2.append(result.rsquared)
print(result.rsquared)
# print(result.pvalues[1:])
plt.title(title + ', r$^2$=%0.2f' % result.rsquared)
