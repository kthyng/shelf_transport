
'''
aggregate connections between significant locations along the coast.
'''

from glob import glob
import numpy as np
import pandas as pd
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

# run in shelf_transport
Files = sorted(glob('calcs/alongcoastconn/conn_in_time/2*.npz'))
# decimal days
t = np.load('calcs/alongcoastconn/conn_in_time/pa_ss/t.npz')['t']

base = 'calcs/alongcoastconn/conn_in_time/between_locs/'
os.makedirs(base, exist_ok=True)
for File in Files:
    # File = Files[0]
    day = File.split('/')[-1][:-4]
    savename = '%s/%s.csv' % (base,day)
    if os.path.exists(savename):
        continue
    d = np.load(File)
    index = [(pd.Timestamp(day) + pd.Timedelta(str(tt) + ' days')).round('H') for tt in t]
    df = pd.DataFrame(index=index)
    mat = d['mat']; d.close()
    # subtract self-connectivity off
    mat[0,:,:] -= np.eye(mat.shape[1])
    for loc0, iboxes0 in locs.items():
        for loc1, iboxes1 in locs.items():
            # don't run if same location
            if loc0 == loc1:
                continue
            df['%s to %s' % (loc0, loc1)] = mat[:,iboxes0,iboxes1[0]:iboxes1[-1]+1].sum(axis=1).sum(axis=1)
    df.to_csv(savename)
