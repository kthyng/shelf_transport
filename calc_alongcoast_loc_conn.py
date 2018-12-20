
'''
aggregate connections between significant locations along the coast.
'''

from glob import glob
import numpy as np
import pandas as pd

# boxes for each location
# bpa = [111,112,113,114,115]; bss = [154,155,156,157,158]

locs = {'Port Mansfield': [79,80, 81, 82, 83],
        'Port Aransas': [111,112,113,114,115],
        "Port O'Connor": [129, 130, 131, 132, 133],
        'Surfside': [154,155,156,157,158],
        'Galveston': [168, 169, 170, 171, 172],
        'Atchafalaya': [237, 238, 239, 240, 241],
        'Terrebonne': [260, 261, 262, 263, 264],
        'Barataria': [274, 275, 276, 277, 278]}

# run in shelf_transport
Files = glob('calcs/alongcoastconn/conn_in_time/2*.npz')
# decimal days
t = np.load('calcs/alongcoastconn/conn_in_time/pa_ss/t.npz')['t']
# dt = t[1] - t[0]  # in days

for File in Files:
    # File = Files[0]
    d = np.load(File)
    day = File.split('/')[-1][:-4]
    index = [(pd.Timestamp(day) + pd.Timedelta(str(tt) + ' days')).round('H') for tt in t]
    df = pd.DataFrame(index=index)
    for loc0, iboxes0 in locs.items():
        for loc1, iboxes1 in locs.items():
            # don't run if same location
            if loc0 == loc1:
                continue
            df['%s to %s' % (loc0, loc1)] = d['mat'][:,iboxes0,iboxes1].sum(axis=1)
