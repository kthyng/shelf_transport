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

# loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
# grid = tracpy.inout.readgrid(loc, usebasemap=True)

# get histogram edges
dstart = np.load('calcs/dispersion/hist/Hstart_bins100.npz')
xe = dstart['xe']; ye = dstart['ye']
dstart.close()

# Where plotting
# plt.figure()
# tracpy.plotting.background(grid)
# plt.plot(xe[45], ye[70], 'o', color='g', ms=15) # shelf eddy region
# plt.plot(xe[15], ye[57], 'o', color='b', ms=15) # Port Aransas bend region

# get time info
dtracks = netCDF.Dataset('tracks/2004-01-01T00gc.nc')
tp = dtracks.variables['tp']
days = (tp-tp[0])/(3600.*24)

# Summer, r=1, just GLAD time: July 20-31, 2012
fnameS = 'calcs/dispersion/hist/D2/r1/D2aveS-GLAD.npz'
if os.path.exists(fnameS):
    D2aveSGLAD = np.load(fnameS)['D2aveS']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/D2/r1/2012-07-[2-3]?T0?_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveSGLAD = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveSGLAD, nnans=nnans)

# Summer, r=5, just GLAD time: July 20-31, 2012
fnameS = 'calcs/dispersion/hist/D2/r5/D2aveS-GLAD.npz'
if os.path.exists(fnameS):
    D2aveS5GLAD = np.load(fnameS)['D2aveS']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/D2/r5/2012-07-[2-3]?T0?_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveS5GLAD = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveS5GLAD, nnans=nnans)

# Summer, r=1
fnameS = 'calcs/dispersion/hist/D2/r1/D2aveS.npz'
if os.path.exists(fnameS):
    D2aveS = np.load(fnameS)['D2aveS']
    daysD2 = np.load(fnameS)['days']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/D2/r1/20??-0[7,8]-??T0?_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveS = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveS)

# Summer, r=5
fnameS = 'calcs/dispersion/hist/D2/r5/D2aveS.npz'
if os.path.exists(fnameS):
    D2aveS5 = np.load(fnameS)['D2aveS']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/D2/r5/20??-0[7,8]-??T0?_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveS5 = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameS, days=days, xe=xe, ye=ye, D2aveS=D2aveS5, nnans=nnans)

# Winter, r=1
fnameW = 'calcs/dispersion/hist/D2/r1/D2aveW.npz'
if os.path.exists(fnameW):
    D2aveW = np.load(fnameW)['D2aveW']
else:
    D2 = np.zeros((99,99,901))
    nnans = 0
    Files = glob.glob('calcs/dispersion/hist/D2/r1/20??-0[1,2]-??T0?_bins100.npz')
    for File in Files:
        print File
        d = np.load(File)
        D2 = np.nansum(np.array((D2, d['D2']*d['nnans'])), axis=0)
        nnans += d['nnans']
        d.close()
    D2aveW = D2/nnans # average over two months of pairs in this histogram bin
    np.savez(fnameW, days=days, xe=xe, ye=ye, D2aveW=D2aveW)



# Plot relative dispersion for LaCasce data vs. different years I have available.
# This is copied from projects/sculp1/dispersion.py
base = '/home/kthyng/projects/sculp1/'
colors = ['.5', '.2', '.4'] # no diff, doturb=1, doturb=2,3
symbols = [':', '-.', '-'] 
lacasce = np.loadtxt(base + 'lacasce_dispersion_Points.txt')
lacasce50 = np.loadtxt(base + 'lacasce_disperion_50days.txt')

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.set_xlabel('Time [days]'); ax.set_ylabel('Mean squared separation distance [km$^2\!$]');
# Skip 2012 because of model output not matching variable set
years = np.array([2003,2004,2005,2006,2007,2008,2009,2010,2011,2013])
d = netCDF.Dataset(base + 'tracks/2003-10-01T00gc.nc')
t = d.variables['tp']
days = (t[0]-t[0,0])/(3600.*24) # time vector in days
d.close()
i25 = find(days<=25)[-1]
D2mean = 0; nnansmean = 0;
for year in years:
    Dname = base + 'calcs/' + str(year) + 'D2.npz'
    D = np.load(Dname)
    D2 = D['D2']; nnans = D['nnans']; D.close()
    D2mean += D2*nnans
    nnansmean += nnans
    ax.semilogy(days[:i25], D2[:i25], color='0.5', linewidth=4, alpha=0.5)  # each year available
ax.semilogy(lacasce[:,0], lacasce[:,1], 'ro', label='LaCasce', markersize=13, alpha=1)  # 25 day lacasce data
# labels for lines
ax.text(19, 2*10**4, 'LaCasce data', color='r', fontsize=18)
ax.text(17, 4*10**2, 'Model comparisons', color='0.5', fontsize=18)
ax.set_ylim(.5, 10**5)
ax.set_xlim(0, 30)
ax.text(0.575, -0.085, 'Data and image from LaCasce and Ohlmann 2003', color='0.2', fontsize=12, transform=ax.transAxes)
fig.savefig('figures/D2/relative_dispersion_comp.pdf', bbox_inches='tight')


## Plot of data and model (lighter, several years) with shelf eddy region added ##
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.set_xlabel('Time [days]'); ax.set_ylabel('Mean squared separation distance [km$^2\!$]');
# Skip 2012 because of model output not matching variable set
years = np.array([2003,2004,2005,2006,2007,2008,2009,2010,2011,2013])
d = netCDF.Dataset(base + 'tracks/2003-10-01T00gc.nc')
t = d.variables['tp']
days = (t[0]-t[0,0])/(3600.*24) # time vector in days
d.close()
i25 = find(days<=25)[-1]
D2mean = 0; nnansmean = 0;
for year in years:
    Dname = base + 'calcs/' + str(year) + 'D2.npz'
    D = np.load(Dname)
    D2 = D['D2']; nnans = D['nnans']; D.close()
    D2mean += D2*nnans
    nnansmean += nnans
    ax.semilogy(days[:i25], D2[:i25], color='0.5', linewidth=4, alpha=0.1)  # each year available
ax.semilogy(lacasce[:,0], lacasce[:,1], 'ro', label='LaCasce', markersize=13, alpha=0.2)  # 25 day lacasce data
# shelf eddy region
ax.semilogy(daysD2, D2aveS[70,45,:], '--', color='g', lw=4) # summer
ax.semilogy(daysD2, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
# labels for lines
ax.text(26, 1.5*10**4, 'Summer', color='g', fontsize=18)
ax.text(26.2, 10**3, 'Winter', color='g', fontsize=18)
ax.set_ylim(.5, 10**5)
ax.set_xlim(0, 30)
ax.text(0.275, -0.085, 'GLAD data from Poje et al, PNAS 2014', color='0.2', fontsize=12, transform=ax.transAxes)
ax.text(0.575, -0.085, 'Data from LaCasce and Ohlmann 2003', color='0.2', fontsize=12, transform=ax.transAxes)
fig.savefig('figures/D2/data-model-shelfregion.pdf', bbox_inches='tight')
fig.savefig('figures/D2/data-model-shelfregion.png', bbox_inches='tight', dpi=300)


## Plot with lacasce and several years of model data. 
## And shelf eddy region winter and summer.
## And GLAD data line.
# load in CARTHE data
# txtc1 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r0pt5.txt', skiprows=1)
txtc2 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r5pt0.txt', skiprows=1)
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.set_xlabel('Time [days]'); ax.set_ylabel('Mean squared separation distance [km$^2\!$]');
# Skip 2012 because of model output not matching variable set
years = np.array([2003,2004,2005,2006,2007,2008,2009,2010,2011,2013])
d = netCDF.Dataset(base + 'tracks/2003-10-01T00gc.nc')
t = d.variables['tp']
days = (t[0]-t[0,0])/(3600.*24) # time vector in days
d.close()
i25 = find(days<=25)[-1]
D2mean = 0; nnansmean = 0;
for year in years:
    Dname = base + 'calcs/' + str(year) + 'D2.npz'
    D = np.load(Dname)
    D2 = D['D2']; nnans = D['nnans']; D.close()
    D2mean += D2*nnans
    nnansmean += nnans
    ax.semilogy(days[:i25], D2[:i25], color='0.5', linewidth=4, alpha=0.1)  # each year available
ax.semilogy(lacasce[:,0], lacasce[:,1], 'ro', label='LaCasce', markersize=13, alpha=0.2)  # 25 day lacasce data
# shelf eddy region
ax.semilogy(daysD2, D2aveS[70,45,:], '--', color='g', lw=4, alpha=0.2) # summer
ax.semilogy(daysD2, D2aveW[70,45,:], '-.', color='g', lw=4, alpha=0.2) # winter
# model to compare with GLAD
# plt.semilogy(daysD2, D2aveSGLAD[70,45,:], '--', color='g', lw=4) # summer, r=1
plt.semilogy(daysD2, D2aveS5GLAD[70,45,:], '--', color='g', lw=4) # summer, r=5
# GLAD
plt.semilogy(txtc2[:,0], txtc2[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# labels for lines
ax.text(26, 3*10**4, 'Summer', color='g', fontsize=16)
ax.text(20, 3*10**4, 'GLAD, r=5km', color='orange', fontsize=16, alpha=0.6)
ax.set_ylim(.5, 10**5)
ax.set_xlim(0, 30)
ax.text(0.575, -0.085, 'Data from LaCasce and Ohlmann 2003', color='0.2', fontsize=12, transform=ax.transAxes)
fig.savefig('figures/D2/data-model-shelfregion-carthe-5km.pdf', bbox_inches='tight')
fig.savefig('figures/D2/data-model-shelfregion-carthe-5km.png', bbox_inches='tight', dpi=300)

# ## Plot all with CARTHE data too and 5km initial separation distance for model ##
# plt.figure(figsize=(12,6))
# # plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
# # plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
# plt.xlabel('Time [days]')
# plt.ylabel('Mean squared separation distance [km$^2\!$]')
# plt.semilogy(days, D2aveSGLAD[70,45,:], '--', color='g', lw=4) # summer, r=1
# plt.semilogy(days, D2aveS5GLAD[70,45,:], '--', color='g', lw=4) # summer, r=5
# # plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
# # plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
# # plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
# plt.semilogy(txtc1[:,0], txtc1[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.semilogy(txtc2[:,0], txtc2[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.savefig('figures/D2/data-model-bothregions-carthe-5km.pdf', bbox_inches='tight')
# plt.savefig('figures/D2/data-model-bothregions-carthe-5km.png', bbox_inches='tight', dpi=300)
# ####



# ## Plot of data and model (lighter) with Port Aransas region added ##
# plt.figure(figsize=(12,6))
# plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
# plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
# plt.xlabel('Time [days]')
# plt.ylabel('Mean squared separation distance [km$^2\!$]')
# # shelf eddy region
# plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
# plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
# plt.savefig('figures/D2/data-model-paregion.pdf', bbox_inches='tight')
# plt.savefig('figures/D2/data-model-paregion.png', bbox_inches='tight', dpi=300)
# ####

# ## Plot of data and model (lighter) with both regions added ##
# plt.figure(figsize=(12,6))
# plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
# plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
# plt.xlabel('Time [days]')
# plt.ylabel('Mean squared separation distance [km$^2\!$]')
# plt.semilogy(days, D2aveS[70,45,:], '--', color='g', lw=4) # summer
# plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
# plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
# plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
# plt.savefig('figures/D2/data-model-bothregions.png', bbox_inches='tight', dpi=300)
# plt.savefig('figures/D2/data-model-bothregions.pdf', bbox_inches='tight')
# ####


# ## Plot all with CARTHE data too ##
# # load in CARTHE data
# txtc1 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r0pt5.txt', skiprows=1)
# txtc2 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r5pt0.txt', skiprows=1)
# plt.figure(figsize=(12,6))
# # plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
# # plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
# plt.xlabel('Time [days]')
# plt.ylabel('Mean squared separation distance [km$^2\!$]')
# plt.semilogy(days, D2aveS[70,45,:], '--', color='g', lw=4) # summer
# # plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
# plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
# # plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
# plt.semilogy(txtc1[:,0], txtc1[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.semilogy(txtc2[:,0], txtc2[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.savefig('figures/D2/data-model-bothregions-carthe.pdf', bbox_inches='tight')
# plt.savefig('figures/D2/data-model-bothregions-carthe.png', bbox_inches='tight', dpi=300)
# ####


# ## Plot all with CARTHE data too and 5km initial separation distance for model ##
# # load in CARTHE data
# txtc1 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r0pt5.txt', skiprows=1)
# txtc2 = np.loadtxt('calcs/dispersion/carthe-fig5-middle-S2-r5pt0.txt', skiprows=1)
# plt.figure(figsize=(12,6))
# # plt.semilogy(txt[:,0], txt[:,1], 'r*', ms=13, zorder=1, alpha=0.4)
# # plt.semilogy(daysold[:i25], D2[:i25], '-', color='0.5', linewidth=4, zorder=1, alpha=0.4)
# plt.xlabel('Time [days]')
# plt.ylabel('Mean squared separation distance [km$^2\!$]')
# plt.semilogy(days, D2aveSGLAD[70,45,:], '--', color='g', lw=4) # summer, r=1
# plt.semilogy(days, D2aveS5GLAD[70,45,:], '--', color='g', lw=4) # summer, r=5
# # plt.semilogy(days, D2aveW[70,45,:], '-.', color='g', lw=4) # winter
# # plt.semilogy(days, D2aveS[57,15,:], '--', color='b', lw=4) # summer
# # plt.semilogy(days, D2aveW[57,15,:], '-.', color='b', lw=4) # winter
# plt.semilogy(txtc1[:,0], txtc1[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.semilogy(txtc2[:,0], txtc2[:,1]**2, '-', color='orange', lw=4, zorder=1, alpha=0.6)
# plt.savefig('figures/D2/data-model-bothregions-carthe-5km.pdf', bbox_inches='tight')
# plt.savefig('figures/D2/data-model-bothregions-carthe-5km.png', bbox_inches='tight', dpi=300)
# ####
