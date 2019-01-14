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



base = 'calcs/alongcoastconn/conn_in_time/'
lnamest = base + 'lines_by_start.npz'
lnameen = base + 'lines_by_end.npz'
# if os.path.exists(lname):
#     d = np.load(lname)
#     startingp = d['startingp']; startingn = d['startingn']
#     endingp = d['endingp']; endingn = d['endingn']; dates = d['dates']; dist = d['dist']
#
# else:
Files = sorted(glob(base + '*.npz'))

t = np.load(Files[0])['t']

startingp = np.empty((len(Files), t.size, 342));
startingn = np.empty((len(Files), t.size, 342))
endingp = np.empty((len(Files), t.size, 342));
endingn = np.empty((len(Files), t.size, 342))
dates = []
i = 0

for i, File in enumerate(Files):

    # print(File)
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
    # i += 1

np.savez(lnamest, t=t, startingp=startingp, startingn=startingn, dates=dates, dist=dist)
np.savez(lnameen, t=t, endingp=endingp, endingn=endingn, dates=dates, dist=dist)
