'''
Plot dispersion in time.
'''

import tracpy
import tracpy.plotting
import netCDF4 as netCDF
import numpy as np
import matplotlib.pyplot as plt
import glob


loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# get histogram edges
dstart = np.load('calcs/dispersion/hist/Hstart_bins100.npz')
xe = dstart['xe']; ye = dstart['ye']
dstart.close()

# Where plotting
plt.figure()
tracpy.plotting.background(grid)
plt.plot(xe[45], ye[70], 'rx')

# get time info
dtracks = netCDF.Dataset('tracks/2004-01-01T00gc.nc')
tp = dtracks.variables['tp']
days = (tp-tp[0])/(3600.*24)

# Summer
D2 = 0
nnans = 0
Files = glob.glob('calcs/dispersion/hist/20??-0[7,8]-??T00_bins100.npz')
for File in Files:
    d = np.load(File)
    # CHANGE BOTH LOOPS TO ACCOUNT FOR NANS
    D2 += d['D2']*d['nnans']
    nnans += d['nnans']
    d.close()
D2aveS = D2/nnans # average over two months of pairs in this histogram bin
np.savez('calcs/dispersion/hist/D2aveS.npz', days=days, xe=xe, ye=ye, D2aveS=D2aveS)

# Winter
D2 = 0
nnans = 0
Files = glob.glob('calcs/dispersion/hist/20??-0[1,2]-??T00_bins100.npz')
for File in Files:
    d = np.load(File)
    D2 += d['D2']*d['nnans']
    nnans += d['nnans']
    d.close()
D2aveW = D2/nnans # average over two months of pairs in this histogram bin

plt.figure()
plt.semilogy(days, D2aveS[70,45,:], 'g', lw=2)
plt.semilogy(days, D2aveW[70,45,:], 'b', lw=2)

# compare with lacasce data
txt = np.loadtxt('/pong/raid/kthyng/projects/horizontal_diffusivity/figures/lacasce_disperion_50days.txt')
plt.semilogy(txt[:,0], txt[:,1], 'r')

# add in model output that matches lacasce data timing for another comparison