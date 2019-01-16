'''
Fix files made in make_conn_plots.run_with_times(), which are called
calcs/alongcoastconn/conn_in_time/2004-01-01T00.npz.
'''

import numpy as np
from glob import glob
import os
import shutil

base = 'calcs/alongcoastconn/conn_in_time'
Files = sorted(glob('%s/2???-??-??T??.npz' % base))

# move these files to another directory to not overwite
newpath = 'calcs/alongcoastconn/conn_in_time/old_BEFOREFIXING/'
os.makedirs(newpath, exist_ok=True)
[shutil.move(File, newpath) for File in Files]

Files = sorted(glob('%s/2???-??-??T??.npz' % newpath))
t = np.load(Files[0])['t']

for File in Files:
    print(File)
    # File = Files[0]
    fname = '%s/%s' % (base,File.split('/')[-1])
    mat = np.load(File)['mat']

    mat[0,range(342),range(342)] = 1
    mat[1:,range(342),range(342)]  = 0

    np.savez(fname, mat=mat, t=t)
