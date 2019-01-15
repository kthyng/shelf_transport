'''
Save along-coast connectivity from single simulation files into
lines.npz files. Copied code from
shelf_transport/notebooks/porta_tofrom_surfside.ipynb, but changing to
not average over months but leave for individual simulations
'''

import numpy as np
import pandas as pd
from glob import glob
import os

# desired resultant resolution of output
# '4hours': preserve time resolution across individual simulation to 4 hourly per each sim
# 'sim': save connectivity by individual simulation
# 'monthly': (this was done in make_conn_plots.run()) save conn by month-year
time_res = '4hours'

vecname = 'calcs/along_coast_wind/coast_vectors.npz'
d = np.load(vecname)
# along-coast distance (km)
dist = d['dist']
d.close()

base = 'calcs/alongcoastconn/conn_in_time/'
os.makedirs(base, exist_ok=True)
lnamest = base + 'lines_by_start_%s.npz' % time_res
lnameen = base + 'lines_by_end_%s.npz' % time_res
Files = sorted(glob(base + '2*.npz'))

t = np.load(Files[0])['t']

if time_res == '4hours':

    startingp = np.empty((len(Files), t.size, 342));
    startingn = np.empty((len(Files), t.size, 342))
    endingp = np.empty((len(Files), t.size, 342));
    endingn = np.empty((len(Files), t.size, 342))

elif time_res == 'sim':

    startingp = np.empty((len(Files), 342));
    startingn = np.empty((len(Files), 342))
    endingp = np.empty((len(Files), 342));
    endingn = np.empty((len(Files), 342))


dates = []

for i, File in enumerate(Files):

    print(File)
    if time_res == '4hours':
        dates.append(pd.Timestamp(File.split('/')[-1][:10]))
        mat = np.load(File)['mat']
        # make downcoast negative
        ix, iy = np.tril_indices(mat.shape[2], k=1)
        mat[:, ix, iy] = -mat[:, ix, iy]
        matp = np.ma.masked_where(mat<0, mat)
        matn = np.ma.masked_where(mat>0, mat)
        startingp[i,:,:] = matp.sum(axis=2)  # gives upcoast alongcoast conn as function of starting position
        endingp[i,:,:] = matp.sum(axis=1)
        startingn[i,:,:] = matn.sum(axis=2)  # gives downcoast alongcoast conn as function of starting position
        endingn[i,:,:] = matn.sum(axis=1)
    elif time_res == 'sim':
        # sums over 30 days of simulation time
        dates.append(pd.Timestamp(File.split('/')[-1][:10]))
        mat = np.load(File)['mat']
        # make downcoast negative
        ix, iy = np.tril_indices(mat.shape[2], k=1)
        mat[:, ix, iy] = -mat[:, ix, iy]
        matp = np.ma.masked_where(mat<0, mat)
        matn = np.ma.masked_where(mat>0, mat)
        # SHOULD THIS BE SUM OR MEAN?
        startingp[i,:] = matp.sum(axis=2).sum(axis=0) # gives upcoast alongcoast conn as function of starting position
        endingp[i,:] = matp.sum(axis=1).sum(axis=0)
        startingn[i,:] = matn.sum(axis=2).sum(axis=0)  # gives downcoast alongcoast conn as function of starting position
        endingn[i,:] = matn.sum(axis=1).sum(axis=0)

np.savez(lnamest, t=t, startingp=startingp, startingn=startingn, dates=dates, dist=dist)
np.savez(lnameen, t=t, endingp=endingp, endingn=endingn, dates=dates, dist=dist)
