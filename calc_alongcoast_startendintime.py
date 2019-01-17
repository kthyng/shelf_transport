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
years = np.arange(2004, 2015)
for year in years:

    # these files were created in make_conn_plots.run_with_times()
    Files = sorted(glob('%s/%s*.npz' % (base, year)))

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

    lnamest = base + 'lines_by_start_%s%i.npz' % (time_res, year)
    lnameen = base + 'lines_by_end_%s%i.npz' % (time_res, year)

    if os.path.exists(lnamest):
        continue

    for i, File in enumerate(Files):

        print(File)
        date = pd.Timestamp(File.split('/')[-1][:10])

        dates.append(date)

        # sums over 30 days of simulation time
        mat = np.load(File)['mat']
        # make downcoast negative
        ix, iy = np.tril_indices(mat.shape[2], k=1)
        #mat is 180 times (in a simulation) x 342 alongcoast boxes x 342 alongcoast boxes
        mat[:, ix, iy] = -mat[:, ix, iy]
        matp = np.ma.masked_where(mat<0, mat)
        matn = np.ma.masked_where(mat>0, mat)

        if time_res == '4hours':
            # subtract or add one at a later accounting point to account for adding 100% connectivity with its own box
            startingp[i,:,:] = matp.sum(axis=2)  # gives upcoast alongcoast conn as function of starting position
            endingp[i,:,:] = matp.sum(axis=1)
            startingn[i,:,:] = matn.sum(axis=2)  # gives downcoast alongcoast conn as function of starting position
            endingn[i,:,:] = matn.sum(axis=1)
        elif time_res == 'sim':
            # SHOULD THIS BE SUM OR MEAN?
            startingp[i,:] = matp.sum(axis=2).sum(axis=0) # gives upcoast alongcoast conn as function of starting position
            endingp[i,:] = matp.sum(axis=1).sum(axis=0)
            startingn[i,:] = matn.sum(axis=2).sum(axis=0)  # gives downcoast alongcoast conn as function of starting position
            endingn[i,:] = matn.sum(axis=1).sum(axis=0)

    np.savez(lnamest, t=t, startingp=startingp, startingn=startingn, dates=dates, dist=dist)
    np.savez(lnameen, t=t, endingp=endingp, endingn=endingn, dates=dates, dist=dist)
