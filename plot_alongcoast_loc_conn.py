
'''
Plot connectivity between locations calculation in calc_alongcoast_loc_conn.py.
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# map to have as subset
merc = cartopy.crs.Mercator()
pc = cartopy.crs.PlateCarree()
land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cartopy.feature.COLORS['land'])
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')


locs = {'Brownsville': [67, 68, 69, 70, 71],
        'Port Mansfield': [79,80, 81, 82, 83],
        'Port Aransas': [111,112,113,114,115],
        "Port O'Connor": [129, 130, 131, 132, 133],
        'Surfside': [154,155,156,157,158],
        'Galveston': [168, 169, 170, 171, 172],
        'Atchafalaya': [237, 238, 239, 240, 241],
        'Terrebonne': [260, 261, 262, 263, 264],
        'Barataria': [274, 275, 276, 277, 278]}

locsorder = ['Barataria', 'Terrebonne', 'Atchafalaya', 'Galveston', 'Surfside',
             "Port O'Connor", 'Port Aransas', 'Port Mansfield', 'Brownsville']

# distance between center box of locations:
# for loc0 in locsorder:
#     for loc1  in locsorder:
#         if loc0 == loc1:
#             continue
#         print('distance from %s to %s is %4.0fkm' % (loc0, loc1, abs(dist[locs[loc0][2]] - dist[locs[loc1][2]])))

colu = '#218983'  # upcoast color
cold = '#cb6863'  # downcoast color

# decimal days
t = np.load('calcs/alongcoastconn/conn_in_time/pa_ss/t.npz')['t']


# read in files
base = 'calcs/alongcoastconn/conn_in_time/between_locs/'
allFiles = sorted(glob('%s/*.csv' % base))
dates = [pd.Timestamp(File.split('/')[-1][:-4]) for File in allFiles]

months = [1,2]
winFiles = [File for date, File in zip(dates, allFiles) if date.month in months]
df = pd.read_csv(winFiles[0], parse_dates=True, index_col=0)

# aggregate files seasonally
fname = base + 'sum_winter.npz'
if not os.path.exists(fname):
    # winFiles = [File for date, File in zip(dates, allFiles) if date.month in months]
    win = np.zeros((len(winFiles), df.columns.size, t.size))  # sim file x locations x 30 days time
    for i, File in enumerate(winFiles):
        df = pd.read_csv(File, parse_dates=True, index_col=0)
        for j, col in enumerate(df.columns):
            win[i,j,:] += df[col]
    np.savez(fname, t=t, dates=dates, win=win, columns=df.columns)
else:
    d = np.load(fname)
    t = d['t']; dates = d['dates']; win = d['win']
    d.close()

months = [7,8]
fname = base + 'sum_summer.npz'
if not os.path.exists(fname):
    sumFiles = [File for date, File in zip(dates, allFiles) if date.month in months]
    sum = np.zeros((len(sumFiles), df.columns.size, t.size))  # sim file x locations x 30 days time
    for i, File in enumerate(sumFiles):
        df = pd.read_csv(File, parse_dates=True, index_col=0)
        for j, col in enumerate(df.columns):
            sum[i,j,:] += df[col]
    np.savez(fname, t=t, dates=dates, sum=sum)
else:
    d = np.load(fname)
    t = d['t']; dates = d['dates']; sum = d['sum']
    d.close()

columns = df.columns

# plot up results: mean and standard deviation for 2 seasons
base = 'figures/alongcoastconn/conn_in_time/between_locs/all_lines'
os.makedirs(base, exist_ok=True)
for loc0 in locs.keys():
    fig, axes = plt.subplots(len(locs)-1, 1, sharex=True, sharey=True, figsize=(8,7))
    fig.suptitle('Transport from %s to ...' % loc0)
    fig.subplots_adjust(top=0.96, left=0.07, right=0.98, hspace=0.35, bottom=0.07)
    axes[-1].set_xlabel('Time [days]')
    axes[0].set_ylim(0, 2)
    axes[0].set_xlim(0, 30)
    i = 0
    for loc1 in locsorder:
        if loc0 == loc1:
            continue
        ax = axes[i]
        col = '%s to %s' % (loc0, loc1)
        icol = np.where(columns == col)[0][0]
        if locs[loc1][0] > locs[loc0][0]:  # loc1 is upcoast from loc0
            color = colu
            # sign = 1
        else:  # loc1 is downcoast from loc0
            color = cold
            # sign = -1
        # winmean = win[:,icol,:].mean(axis=0)
        # winstd = win[:,icol,:].std(axis=0)
        # ax.fill_between(t, winmean-winstd, winmean+winstd, color=color, alpha=0.3)
        ax.plot(t, win[:,icol,:].T, color=color, alpha=0.2, lw=2);
        ax.plot(t, sum[:,icol,:].T, color=color, alpha=0.6, lw=1);
        ax.text(0.8, 0.8, loc1, transform=ax.transAxes)
        if i == 0:
            axes[i].text(0.01, 0.7, 'winter', color=color, alpha=0.4, transform=axes[i].transAxes)
            axes[i].text(0.15, 0.7, 'summer', color=color, alpha=0.8, transform=axes[i].transAxes)
        if i == 1:
            axes[i].text(0.01, 0.7, 'downcoast', color=cold, alpha=0.6, transform=axes[i].transAxes)
            axes[i].text(0.15, 0.7, 'upcoast', color=colu, alpha=0.6, transform=axes[i].transAxes)
        i += 1
    fig.savefig('%s/%s.png' % (base,loc0), bbox_inches='tight')
    fig.savefig('%s/%s.pdf' % (base,loc0), bbox_inches='tight')
    plt.close(fig)

# calculate new arrays where if any connectivity has occurred, it is marked as
# yes for the rest of the 30 days
winconn = np.zeros_like(win)
for icol, col in enumerate(columns):
    for irow, row in enumerate(win[:,icol,:]):  # loops through Files
        try:
            i1st = np.where(row > 0)[0][0]
            winconn[irow,icol,i1st:] = 1
        except:
            pass
sumconn = np.zeros_like(sum)
for icol, col in enumerate(columns):
    for irow, row in enumerate(sum[:,icol,:]):  # loops through Files
        try:
            i1st = np.where(row > 0)[0][0]
            sumconn[irow,icol,i1st:] = 1
        except:
            pass

vecname = 'calcs/along_coast_wind/coast_vectors.npz'
dist = np.load(vecname)['dist']

# save plot properties
params = {}
for loc0 in locs.keys():
    params[loc0] = {}
params['Port Mansfield']['axes'] = [0.15, 0.35, 0.3, 0.8]
params['Terrebonne']['axes'] = [0.15, 0.1, 0.3, 0.8]
params['Barataria']['axes'] = [0.15, -0.1, 0.3, 0.8]
params['Atchafalaya']['axes'] = [0.15, 0.01, 0.3, 0.8]
params['Galveston']['axes'] = [0.12, -0.1, 0.3, 0.8]
params['Surfside']['axes'] = [0.105, -0.225, 0.27, 0.8]
params["Port O'Connor"]['axes'] = [0.105, 0.4275, 0.3, 0.8]
params['Port Aransas']['axes'] = [0.15, 0.4, 0.3, 0.8]
params['Brownsville']['axes'] = [0.105, 0.35, 0.3, 0.8]

# axes index and left (x) starting location for winter text (summer to right)
params['Port Mansfield']['text'] = [2, 0.02]
params['Terrebonne']['text'] = [2, 0.02]
params['Barataria']['text'] = [2, 0.02]
params['Atchafalaya']['text'] = [0, 0.02]
params['Galveston']['text'] = [2, 0.02]
params['Surfside']['text'] = [2, 0.02]
params["Port O'Connor"]['text'] = [4, 0.02]
params['Port Aransas']['text'] = [3, 0.02]
params['Brownsville']['text'] = [0, 0.51]

base = 'figures/alongcoastconn/conn_in_time/between_locs/sum_lines'
os.makedirs(base, exist_ok=True)
# save no-transport times for table
df = pd.DataFrame(index=np.arange(0,9*8), columns=['from', 'to', 'dist [km]', 'min time winter [days]',
                           'max speed winter [m/s]', 'min time summer [days]',
                           'max speed summer [m/s]'])
j = 0
for loc0 in locsorder:
    # run through bay options once to see how many there are and store time series
    res = {}  # save results in dictionary
    for loc1 in locsorder:
        if loc0 == loc1:
            continue
        col = '%s to %s' % (loc0, loc1)
        icol = np.where(columns == col)[0][0]
        winmean = ((winconn[:,icol,:]>0).sum(axis=0)/(win.shape[0]))
        winmean[winmean==0] = np.nan
        winstd = ((winconn[:,icol,:]>0).std(axis=0))
        winstdm = (winmean-winstd)*100; winstdp = (winmean+winstd)*100
        winstdm[winstdm<=0] = 0; winstdp[winstdp>100] = 100
        summean = ((sumconn[:,icol,:]>0).sum(axis=0)/(sum.shape[0]))
        summean[summean==0] = np.nan
        sumstd = ((sumconn[:,icol,:]>0).std(axis=0))
        sumstdm = (summean-sumstd)*100; sumstdp = (summean+sumstd)*100
        sumstdm[sumstdm<=0] = 0; sumstdp[sumstdp>100] = 100

        # # write no-transport times to file
        # inotransport = np.where(np.isnan(winmean))[0][-1]
        # f.write('From %s to %s in winter: %2.1f days\n' % (loc0, loc1, t[inotransport]))
        # inotransport = np.where(np.isnan(summean))[0][-1]
        # f.write('From %s to %s in summer: %2.1f days\n' % (loc0, loc1, t[inotransport]))

        # also save in dataframe
        df['from'].iloc[j] = loc0
        df['to'].iloc[j] = loc1
        # distance from center box to center box
        df['dist [km]'].iloc[j] = abs(dist[locs[loc0][2]] - dist[locs[loc1][2]])
        inotransport = np.where(np.isnan(winmean))[0][-1]
        if t[inotransport] == 30:
            tt = np.nan
        else:
            tt = t[inotransport]
        df['min time winter [days]'].iloc[j] = tt
        df['max speed winter [m/s]'].iloc[j] = df['dist [km]'].iloc[j]*1000/(tt*86400)
        inotransport = np.where(np.isnan(summean))[0][-1]
        if t[inotransport] == 30:
            tt = np.nan
        else:
            tt = t[inotransport]
        df['min time summer [days]'].iloc[j] = tt
        df['max speed summer [m/s]'].iloc[j] = df['dist [km]'].iloc[j]*1000/(tt*86400)
        j += 1

        # don't use if all mean values are nan
        if (np.isnan(winmean).sum() == winmean.size) and \
            (np.isnan(summean).sum() == summean.size):
            pass
        else:  # save
            res[loc1] = {}
            res[loc1]['winmean'] = winmean
            res[loc1]['summean'] = summean
            res[loc1]['winstdm'] = winstdm
            res[loc1]['sumstdm'] = sumstdm
            res[loc1]['winstdp'] = winstdp
            res[loc1]['sumstdp'] = sumstdp

    nrows = len(res.keys())
    fig, axes = plt.subplots(nrows, 1, sharex=True, sharey=True, figsize=(5,nrows*0.85))
    fig.subplots_adjust(top=0.95, left=0.1, right=0.98, hspace=0.35, bottom=0.07)
    axes[-1].set_xlabel('Time [days]')
    # axes[-1].set_ylabel('Likelihood any connectivity [%]', horizontalalignment='left')
    irow = int(np.floor(nrows/2))
    axes[irow].text(-0.125, 0.5, 'Likelihood any connectivity [%]',
                    verticalalignment='center', rotation=90,
                    transform=axes[irow].transAxes)
    axes[0].text(0.01, 1.1, 'Transport from %s to ...' % loc0, transform=axes[0].transAxes, fontsize=12)
    axes[0].set_ylim(0, 100)
    axes[0].set_xlim(0, 30)
    i = 0
    for loc1 in locsorder:
        if loc0 == loc1:
            continue
        # only use to-location if has non-nan results
        if not loc1 in res.keys():
            continue
        ax = axes[i]
        if locs[loc1][0] > locs[loc0][0]:  # loc1 is upcoast from loc0
            color = colu
        else:  # loc1 is downcoast from loc0
            color = cold
        ax.plot(t, res[loc1]['winmean']*100, ls='-', alpha=0.6, lw=3, color=color)
        ax.fill_between(t, res[loc1]['winstdm'], res[loc1]['winstdp'], alpha=0.15, color=color)
        ax.plot(t, res[loc1]['summean']*100, ls='--', alpha=0.6, lw=3, color=color)
        ax.fill_between(t, res[loc1]['sumstdm'], res[loc1]['sumstdp'], alpha=0.15, color=color, hatch='/')
        ax.text(0.775, 1.05, loc1, transform=ax.transAxes)
        ax.grid(True)
        if i == params[loc0]['text'][0]:
            axes[i].text(params[loc0]['text'][1], 0.75, 'winter $-$', color=color, alpha=0.8, transform=axes[i].transAxes)
            axes[i].text(params[loc0]['text'][1]+0.16, 0.75, 'summer $--$', color=color, alpha=0.8, transform=axes[i].transAxes)
        i += 1

    # add map
    # fig = plt.figure(figsize=(7, 4))
    axm = fig.add_axes(params[loc0]['axes'], projection=merc)
    axm.set_frame_on(False) # kind of like it without the box
    axm.set_extent([-98, -88.5, 25.5, 30.3], pc)
    gl = axm.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_left = False
    axm.add_feature(land_10m, facecolor='0.85')
    # axm.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
    axm.add_feature(cartopy.feature.BORDERS, facecolor='0.8')
    axm.add_feature(states_provinces, edgecolor='gray')
    # add location
    pts_u_boxes = np.load('calcs/coastpaths_pts.npz')['pts_u_boxes']
    pt = pts_u_boxes[locs[loc0][2]][0]
    axm.plot(*pt, '*', color='yellow', markersize=10, transform=pc, mec='k')
    # connectivity locations
    for loc1 in locsorder:
        if loc0 == loc1:
            continue
        # plot as empty if has nan results
        if not loc1 in res.keys():
            pts_u_boxes = np.load('calcs/coastpaths_pts.npz')['pts_u_boxes']
            pt = pts_u_boxes[locs[loc1][2]][0]
            axm.plot(*pt, 'o', color='none', markersize=5, transform=pc, mec='k')
            continue
        # ax = axes[i]
        # col = '%s to %s' % (loc0, loc1)
        # icol = np.where(columns == col)[0][0]
        if locs[loc1][0] > locs[loc0][0]:  # loc1 is upcoast from loc0
            color = colu
            # sign = 1
        else:  # loc1 is downcoast from loc0
            color = cold
        pts_u_boxes = np.load('calcs/coastpaths_pts.npz')['pts_u_boxes']
        pt = pts_u_boxes[locs[loc1][2]][0]
        axm.plot(*pt, 'o', color=color, markersize=5, transform=pc)#, mec='k')

    fig.savefig('%s/%s.png' % (base,loc0), bbox_inches='tight')
    fig.savefig('%s/%s.pdf' % (base,loc0), bbox_inches='tight')
    plt.close(fig)
# f.close()
df.to_csv('calcs/alongcoastconn/conn_in_time/between_locs/nums.csv')

#
# fig = plt.figure(figsize=(7, 4))
# ax = fig.add_axes([0.1, 0.01, 0.87, 0.95], projection=merc)
# ax.set_frame_on(False) # kind of like it without the box
# ax.set_extent([-98, -88.5, 25.5, 30.3], pc)
# gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
# # the following two make the labels look like lat/lon format
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.xlabels_bottom = False  # turn off labels where you don't want them
# gl.ylabels_right = False
# ax.add_feature(land_10m, facecolor='0.8')
# ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
# ax.add_feature(cartopy.feature.BORDERS, facecolor='0.8')
# ax.add_feature(states_provinces, edgecolor='gray')
#
# # add location
# pts_u_boxes = np.load('calcs/coastpaths_pts.npz')['pts_u_boxes']
# pt = pts_u_boxes[locs['Terrebonne'][2]][0]
# ax.plot(*pt, '*', color='yellow', markersize=20, transform=pc, mec='k')



# colored matrix of speeds
df = pd.read_csv('calcs/alongcoastconn/conn_in_time/between_locs/nums.csv')
import seaborn as sns
import cmocean.cm as cmo
speedswin = df.pivot("from", "to", "max speed winter [m/s]")
speedswin = speedswin.loc[locsorder,locsorder[::-1]]
speedssum = df.pivot("from", "to", "max speed summer [m/s]")
speedssum = speedssum.loc[locsorder,locsorder[::-1]]
vmin = min((speedswin.min().min(), speedssum.min().min()))
vmax = max((speedswin.max().max(), speedssum.max().max()))
# gray diagonal
dd = 0.75
I = np.eye(len(locsorder))[::-1]*dd
mask = dd - I

f, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
cbar_ax = f.add_axes([.91, .355, .015, .6])
i=0
sns.heatmap(I, mask=mask, ax=ax[0], cmap=cmo.gray, vmin=0, vmax=1, square=False, cbar=False, xticklabels=False, yticklabels=False)  # gray diagonal
sns.heatmap(speedswin, annot=True, linewidths=.5, ax=ax[0], cmap=cmo.speed, vmin=vmin, vmax=vmax, square=False, cbar_kws={"label": "max speed [m/s]"}, cbar=i == 0, cbar_ax=None if i else cbar_ax)
ax[0].set_title('Winter')
ax[0].set_ylabel('Starting location')
ax[0].set_xlabel('Ending location')
ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation = 45, ha='right')#, fontsize = 8)
i=1
sns.heatmap(I, mask=mask, ax=ax[1], cmap=cmo.gray, vmin=0, vmax=1, square=False, cbar=False, xticklabels=False, yticklabels=False)  # gray diagonal
sns.heatmap(speedssum, annot=True, linewidths=.5, ax=ax[1], cmap=cmo.speed, vmin=vmin, vmax=vmax, square=False, cbar_kws={"label": "max speed [m/s]"}, cbar=i == 0, cbar_ax=None if i else cbar_ax)
ax[1].set_title('Summer')
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 45, ha='right')#, fontsize = 8)
ax[1].set_ylabel('')
ax[1].set_xlabel('Ending location')
f.tight_layout(rect=[0, 0, .9, 1])
f.savefig('figures/alongcoastconn/conn_in_time/between_locs/speeds.png', bbox_inches='tight')
f.savefig('figures/alongcoastconn/conn_in_time/between_locs/speeds.pdf', bbox_inches='tight')


# colored matrix of times (smaller and prettier than table)
timeswin = df.pivot("from", "to", "min time winter [days]")
timeswin = timeswin.loc[locsorder,locsorder[::-1]]
timessum = df.pivot("from", "to", "min time summer [days]")
timessum = timessum.loc[locsorder,locsorder[::-1]]
vmin = min((timeswin.min().min(), timessum.min().min()))
vmax = max((timeswin.max().max(), timessum.max().max()))
# gray diagonal
dd = 0.75
I = np.eye(len(locsorder))[::-1]*dd
mask = dd - I

f, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
cbar_ax = f.add_axes([.91, .355, .015, .6])
i=0
sns.heatmap(I, mask=mask, ax=ax[0], cmap=cmo.gray, vmin=0, vmax=1, square=False, cbar=False, xticklabels=False, yticklabels=False)  # gray diagonal
sns.heatmap(timeswin, annot=True, linewidths=.5, ax=ax[0], cmap=cmo.tempo_r, vmin=vmin, vmax=vmax, square=False, cbar_kws={"label": "min travel time [days]"}, cbar=i == 0, cbar_ax=None if i else cbar_ax)
ax[0].set_title('Winter')
ax[0].set_ylabel('Starting location')
ax[0].set_xlabel('Ending location')
ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation = 45, ha='right')#, fontsize = 8)
i=1
sns.heatmap(I, mask=mask, ax=ax[1], cmap=cmo.gray, vmin=0, vmax=1, square=False, cbar=False, xticklabels=False, yticklabels=False)  # gray diagonal
sns.heatmap(timessum, annot=True, linewidths=.5, ax=ax[1], cmap=cmo.tempo_r, vmin=vmin, vmax=vmax, square=False, cbar_kws={"label": "min travel time [days]"}, cbar=i == 0, cbar_ax=None if i else cbar_ax)
ax[1].set_title('Summer')
ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation = 45, ha='right')#, fontsize = 8)
ax[1].set_ylabel('')
ax[1].set_xlabel('Ending location')
f.tight_layout(rect=[0, 0, .9, 1])
f.savefig('figures/alongcoastconn/conn_in_time/between_locs/times.png', bbox_inches='tight')
f.savefig('figures/alongcoastconn/conn_in_time/between_locs/times.pdf', bbox_inches='tight')
