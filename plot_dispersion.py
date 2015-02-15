'''
Plot dispersion in time.
'''

import tracpy
import tracpy.plotting
import netCDF4 as netCDF
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pdb
from matplotlib.mlab import find
import matplotlib as mpl

mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# get histogram edges
dstart = np.load('calcs/dispersion/hist/Hstart_bins100.npz')
xe = dstart['xe']; ye = dstart['ye']
dstart.close()

# Where plotting
plt.figure()
tracpy.plotting.background(grid)
plt.plot(xe[45], ye[70], 'o', color='g', ms=15) # shelf eddy region
plt.plot(xe[15], ye[57], 'o', color='b', ms=15) # Port Aransas bend region

# get time info
dtracks = netCDF.Dataset('tracks/2004-01-01T00gc.nc')
tp = dtracks.variables['tp']
days = (tp-tp[0])/(3600.*24)

# Summer
fnameS = 'calcs/dispersion/hist/D2aveS.npz'
if os.path.exists(fnameS):
    D2aveS = np.load(fnameS)['D2aveS']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/20??-0[7,8]-??T00_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveS = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveS)

# Winter
fnameW = 'calcs/dispersion/hist/D2aveW.npz'
if os.path.exists(fnameW):
    D2aveW = np.load(fnameW)['D2aveW']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/20??-0[1,2]-??T00_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveW = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveW)


## Plot of LaCasce data and model output ##
plt.figure(figsize=(12,6))
# compare with lacasce data, 25 days data
txt = np.loadtxt('/pong/raid/kthyng/projects/shelf_transport/lacasce_dispersion_Points.txt')
plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1)
# add in model output that matches lacasce data timing for another comparison
base = '/pong/raid/kthyng/projects/horizontal_diffusivity/tracks/'
run = base + 'doturb0_ah0'
Dname = os.path.join(run, 'D2overall.npz')
D = np.load(Dname)
D2 = D['D2']; t = D['t']; D.close()
daysold = (t-t[0])/(3600.*24) # time vector in days
# doturb = int(run.split('doturb')[1][0]) # this gives the doturb value for plotting
# if doturb == 3: doturb=2;
i25 = find(daysold<=25)[-1]
plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1)
plt.xlim(0,30)
plt.xlabel('Time [days]')
plt.ylabel('Mean squared separation distance [km$^2\!$]')
plt.savefig('figures/D2/data-model.pdf', bbox_inches='tight')
####

## Plot of data and model (lighter) with shelf eddy region added ##
plt.figure(figsize=(12,6))
plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
plt.xlabel('Time [days]')
plt.ylabel('Mean squared separation distance [km$^2\!$]')
# shelf eddy region
plt.semilogy(days, D2aveS[70,45,:], '--', color='g', lw=4) # summer
plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
plt.savefig('figures/D2/data-model-shelfregion.pdf', bbox_inches='tight')
####

## Plot of data and model (lighter) with Port Aransas region added ##
plt.figure(figsize=(12,6))
plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
plt.xlabel('Time [days]')
plt.ylabel('Mean squared separation distance [km$^2\!$]')
# shelf eddy region
plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
plt.savefig('figures/D2/data-model-paregion.pdf', bbox_inches='tight')
####

## Plot of data and model (lighter) with both regions added ##
plt.figure(figsize=(12,6))
plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
plt.xlabel('Time [days]')
plt.ylabel('Mean squared separation distance [km$^2\!$]')
plt.semilogy(days, D2aveS[70,45,:], '--', color='g', lw=4) # summer
plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
plt.savefig('figures/D2/data-model-bothregions.pdf', bbox_inches='tight')
####


## Plot all with CARTHE data too ##
# load in CARTHE data
txtc1 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r0pt5.txt', skiprows=1)
txtc2 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r5pt0.txt', skiprows=1)
plt.figure(figsize=(12,6))
plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
plt.xlabel('Time [days]')
plt.ylabel('Mean squared separation distance [km$^2\!$]')
plt.semilogy(days, D2aveS[70,45,:], '--', color='g', lw=4) # summer
plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
plt.semilogy(txtc1[:,0], txtc1[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
plt.semilogy(txtc2[:,0], txtc2[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
plt.savefig('figures/D2/data-model-bothregions-carthe.pdf', bbox_inches='tight')
####
