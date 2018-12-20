
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
# CURRENTLY WINTER ONLY
for loc0 in locs.keys():
    fig, axes = plt.subplots(len(locs)-1, 1, sharex=True, sharey=True, figsize=(10,8))
    fig.suptitle('Transport from %s to ...' % loc0)
    axes[-1].set_xlabel('Time [days]')
    i = 0
    for loc1 in locs.keys():
        if loc0 == loc1:
            continue
        ax = axes[i]
        col = '%s to %s' % (loc0, loc1)
        icol = np.where(columns == col)[0][0]
        if locs[loc1][0] > locs[loc0][0]:  # loc1 is upcoast from loc0
            color = colu
        else:  # loc1 is downcoast from loc0
            color = cold
        ax.plot(t, win[:,icol,:].mean(axis=0), color=color)
        ax.text(0.8, 0.9, loc1, transform=ax.transAxes)
        i += 1
