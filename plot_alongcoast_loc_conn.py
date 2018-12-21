
'''
Plot connectivity between locations calculation in calc_alongcoast_loc_conn.py.
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os


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

colu = '#218983'  # upcoast color
cold = '#cb6863'  # downcoast color

# decimal days
t = np.load('calcs/alongcoastconn/conn_in_time/pa_ss/t.npz')['t']


# read in files
base = 'calcs/alongcoastconn/conn_in_time/between_locs/'
allFiles = sorted(glob('%s/*.csv' % base))
dates = [pd.Timestamp(File.split('/')[-1][:-4]) for File in allFiles]

# aggregate files seasonally
months = [1,2]
fname = base + 'sum_winter.npz'
if not os.path.exists(fname):
    winFiles = [File for date, File in zip(dates, allFiles) if date.month in months]
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


# plot up results: mean and standard deviation for 2 seasons
for loc0 in locs.keys():
    fig, axes = plt.subplots(len(locs)-1, 1, sharex=True, sharey=True, figsize=(8,7))
    fig.suptitle('Transport from %s to ...' % loc0)
    fig.subplots_adjust(top=0.96, left=0.07, right=0.98)
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
        i += 1


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


for loc0 in locs.keys():
    fig, axes = plt.subplots(len(locs)-1, 1, sharex=True, sharey=True, figsize=(8,7))
    fig.suptitle('Transport from %s to ...' % loc0)
    fig.subplots_adjust(top=0.95, left=0.07, right=0.98, bottom=0.05)
    axes[-1].set_xlabel('Time [days]')
    axes[0].set_ylim(0, 50)
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
        ax.plot(t, ((winconn[:,icol,:]>0).sum(axis=0).cumsum()/(win.shape[0]*win.shape[2]))*100, ls='-', alpha=0.4, lw=3, color=color)
        ax.plot(t, ((sumconn[:,icol,:]>0).sum(axis=0).cumsum()/(sum.shape[0]*sum.shape[2]))*100, ls='--', lw=3, alpha=0.4, color=color)
        ax.text(0.025, 0.8, loc1, transform=ax.transAxes)
        ax.grid(True)
        i += 1
