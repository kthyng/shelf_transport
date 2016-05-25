'''
Plot a set of drifters in map view, in time.
'''

import matplotlib.pyplot as plt
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
from matplotlib.mlab import find
import pdb
import numpy as np
import matplotlib as mpl
import os
import tracpy.calcs
import glob
import bisect
from datetime import datetime, timedelta
import op
import matplotlib.patches as patches
import cmocean as cm

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

# read in grid
# loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# loc = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# loc = '/home/kthyng/shelf/grid.nc'
loc = '/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_*.nc'
grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
# grid = tracpy.inout.readgrid(grid_filename, usebasemap=True, llcrnrlat=22.85, llcrnrlon=-97.9, urcrnrlat=30.5)
proj = tracpy.tools.make_proj('nwgom', usebasemap=True, llcrnrlat=27, llcrnrlon=-96, urcrnrlon=-90)
grid = tracpy.inout.readgrid(grid_filename, proj)

# whether to do tails on drifters or not (don't with low decimation)
dotails = False  # True or False
dostreaklines = False  # for when not doing old tails
dowind = False  # whether to plot wind or not
docurrents = False  # whether or not to plot the surface currents
doriver = False  # whether or not to plot river
dopong = False  # whether to label with pong.tamu.edu
docontour = False  # whether or not to overlay an isobath

# Read in drifter tracks
dd = 1 # 500 # drifter decimation
startdate = '2008-07-14T00'  # 11-14 end in 4
# startdate = '2005-02-14T00'  # 11-14 end in 4

# m = netCDF.Dataset(loc)
m = netCDF.MFDataset(loc)
year = int(startdate.split('-')[0])
month = int(startdate.split('-')[1])
day = int(startdate.split('-')[2].split('T')[0])
hour = int(startdate.split('-')[2].split('T')[1])

## Wind forcing ##
if dowind:
    # There are multiple file locations
    if year <= 2012:
        w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
    elif year == 2013:
        w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2013.nc')
    elif year == 2014:
        w = netCDF.Dataset('/rho/raid/home/kthyng/txla/txla_wind_narr_2014.nc')

    # w = netCDF.Dataset('/atch/raid1/zhangxq/Projects/narr_txla/txla_blk_narr_' + str(year) + '.nc')
    # Wind time period to use
    unitsWind = (w.variables['time'].units).replace('/','-')
    datesWind = netCDF.num2date(w.variables['time'][:], unitsWind)
    # datesWind = datesModel
    wdx = 30; wdy = 50 # in indices
    ##

    # to rotate wind vectors
    anglev = m.variables['angle'][:]

    def rot2d(x, y, ang):
        '''rotate vectors by geometric angle'''
        xr = x*np.cos(ang) - y*np.sin(ang)
        yr = x*np.sin(ang) + y*np.cos(ang)
        return xr, yr

xr = np.asanyarray(grid.x_rho.T, order='C')
yr = np.asanyarray(grid.y_rho.T, order='C')
xpsi = np.asanyarray(grid.x_psi.T, order='C')
ypsi = np.asanyarray(grid.y_psi.T, order='C')


# Model time period to use
units = m.variables['ocean_time'].units
starttime = netCDF.date2num(datetime(year, 1, 1, 4, 0, 0), units)
if year==2014:
    endtime = netCDF.date2num(datetime(year, 9, 30, 20, 0, 0), units)
else:
    endtime = netCDF.date2num(datetime(year+1, 1, 1, 4, 0, 0), units)
if docontour:
    dt = m.variables['ocean_time'][1] - m.variables['ocean_time'][0] # 4 hours in seconds
    ts = np.arange(starttime, endtime, dt)
    itshift = find(starttime==m.variables['ocean_time'][:]) # shift to get to the right place in model output
    datesModel = netCDF.num2date(m.variables['ocean_time'][:], units)
    plotdates = netCDF.num2date(ts, units)
# current arrows
cdx = 7; cdy = 11  # in indices

if year == 2014:
    monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,10)]
else:
    monthdates = [datetime(year, month, 1, 0, 0, 0) for month in np.arange(1,13)]

if doriver:
    ## River forcing ##
    r1 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012.nc') # use for through 2011
    r2 = netCDF.Dataset('/rho/raid/home/kthyng/txla/TXLA_river_4dyes_2012_2014.nc') # use for 2012-2014
    # River timing
    tr1 = r1.variables['river_time']
    tunitsr1 = tr1.units
    # interpolate times for this data file since at the 12 hours mark instead of beginning of the day
    tr1 = op.resize(tr1, 0)
    datesr1 = netCDF.num2date(tr1[:], tunitsr1)
    tr2 = r2.variables['river_time']
    datesr2 = netCDF.num2date(tr2[:], tr2.units)
    # all of river input
    Q1 = np.abs(r1.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
    # interpolate this like for time
    Q1 = op.resize(Q1, 0)
    Q2 = np.abs(r2.variables['river_transport'][:]).sum(axis=1)*2.0/3.0
    # Combine river info into one dataset
    iend1 = find(datesr1<datetime(2012,1,1,0,0,0))[-1] # ending index for file 1
    tRiver = np.concatenate((tr1[:iend1], tr2[:]), axis=0)
    datesRiver = np.concatenate((datesr1[:iend1], datesr2))
    R = np.concatenate((Q1[:iend1], Q2))
    r1.close(); r2.close()
    # start and end indices in time for river discharge
    itstartRiver = bisect.bisect_left(datesRiver, datetime(year, 1, 1, 0, 0, 0))
    if year == 2014:
        itendRiver = bisect.bisect_left(datesRiver, datetime(year, 9, 30, 20, 0, 0))
    else:
        itendRiver = bisect.bisect_left(datesRiver, datetime(year+1, 1, 1, 0, 0, 0))
    # ticks for months on river discharge
    mticks = [bisect.bisect_left(datesRiver, monthdate) for monthdate in np.asarray(monthdates)]
    mticknames = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    ##


d = netCDF.Dataset('tracks/' + startdate + 'gc.nc')
xg = d.variables['xg'][::dd,:]
yg = d.variables['yg'][::dd,:]
ind = (xg == -1)
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
xp[ind] = np.nan; yp[ind] = np.nan
if d.variables['tp'].ndim == 1:
    tp = d.variables['tp'][:]
elif d.variables['tp'].ndim > 1:
    # index of a drifter that lasts the whole time
    itp = np.argmax(np.max(d.variables['tp'][:], axis=1))
    tp = d.variables['tp'][itp,:]

units = d.variables['tp'].units
# d.close()

# txla output
## Model output ##
# year = int(startdate.split('-')[0])
# if year <= 2013:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_????.nc'))
# elif year == 2014:
#     currents_filenames = np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_??.nc'))
# nc = netCDF.MFDataset(currents_filenames)
if docontour:
    currents_filenames = loc
    nc = netCDF.Dataset(currents_filenames)
    datestxla = netCDF.num2date(nc.variables['ocean_time'][:], nc.variables['ocean_time'].units)

if not dotails:
    
    # Find indices of drifters by starting depth
    depthp = tracpy.calcs.Var(xg[:,0], yg[:,0], tp, 'h', netCDF.Dataset(grid_filename)) # starting depths of drifters
    # near-shore
    # ind10 = depthp<=10
    ind20 = depthp<=20
    ind50 = (depthp>20)*(depthp<=50)
    ind100 = (depthp>50)*(depthp<=100)
    # offshore
    ind500 = (depthp>100)*(depthp<=500)
    ind3500 = depthp>500

    # colors for drifters
    rgb = cm.cm.temp_r(np.linspace(0, 1, 6))[:-1, :]
    # rgb = plt.cm.get_cmap('winter_r')(np.linspace(0,1,6))[:-1,:3] # skip last entry where it levels off in lightness

    # to plot colorbar
    cmapnew = cm.tools.cmap(rgb)  # cuts off dark end piece since interferes with isobath and isohaline
    gradient = np.linspace(0, 1, 6)[:-1]
    gradient = np.vstack((gradient, gradient))

    ms = 1.5 # markersize


# Plot drifters, starting 5 days into simulation
# 2 days for one tail, 3 days for other tail
# t = tp-tp[0]
days = (tp-tp[0])/(3600.*24)
dates = netCDF.num2date(tp, units)
# Find indices relative to present time
i5daysago = 0 # keeps track of index 5 days ago
i2daysago = find(days>=3)[0] # index for 2 days ago, which starts as 3 days in
if dotails:
    i5days = find(days>=5)[0] # index for 5 days in
else:
    i5days = 0 # start at the beginning
nt = tp.size # total number of time indices
# for i in np.arange(0,nt+1,5):

dirname = 'figures/drifters/dd' + str(dd) + '/' + startdate + 'smalldomain'
if not os.path.exists(dirname):
    os.makedirs(dirname)

for i in np.arange(i5days,nt+1,5):

    fname = dirname + '/' + dates[i].isoformat()[:-6] + '.png'

    if docontour:
        itxla = np.where(datestxla==dates[i])[0][0] # find time index to use for model output
        salt = nc.variables['salt'][itxla,-1,:,:] # surface salinity

    if os.path.exists(fname):
        # Update indices
        i5daysago += 5
        i2daysago += 5
        continue

    # print fname

    # pdb.set_trace()
    if dowind:
        itwind = bisect.bisect_left(datesWind, dates[i]) # index for wind at this time
    if docontour:
        itmodel = bisect.bisect_left(datesModel, dates[i]) # index for model output at this time
    if doriver:
        itriver = bisect.bisect_left(datesRiver, dates[i]) # index for river at this time

    # Plot background
    fig = plt.figure(figsize=(10.1, 8.4), dpi=300)  # dpi=150)
    ax = fig.add_axes([0.06, 0.00, 0.93, 0.97])
    ax.set_frame_on(False) # kind of like it without the box
    tracpy.plotting.background(grid=grid, ax=ax, outline=[1,1,0,1], mers=np.arange(-97, -87), merslabels=[0, 0, 1, 0], pars=np.arange(23, 32), hlevs=[100], col='0.2')

    # # Label isobaths
    # ax.text(0.85, 0.865, '10 m', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    # ax.text(0.88, 0.862, '20', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    # ax.text(0.87, 0.835, '50', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    # ax.text(0.89, 0.825, '100', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)
    # ax.text(0.9, 0.803, '450', transform=ax.transAxes, fontsize=9, color='0.4', rotation=45)

    # Date
    date = dates[i].strftime('%Y %b %02d %H:%M')
    # date = datesModel[itmodel].strftime('%Y %b %02d %H:%M')
    ax.text(0.35, 0.425, date, fontsize=18, color='0.2', transform=ax.transAxes, 
                bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))

    if dopong:
        # PONG
        ax.text(0.6, 0.97, 'pong.tamu.edu', fontsize=12, transform=ax.transAxes, color='0.3')

    if dotails:
        # Plot 5 days ago to 2 days ago
        ax.plot(xp[:,i5daysago:i2daysago].T, yp[:,i5daysago:i2daysago].T, color='0.6', lw=2)

        # Plot 0-2 day tail
        ax.plot(xp[:,i2daysago:i].T, yp[:,i2daysago:i].T, color='0.3', lw=3)

        # Plot drifter locations
        ax.plot(xp[:,i].T, yp[:,i].T, 'o', color='r', ms=10)

    else:

        if dostreaklines:

            istart = i-(24./3*4*0.25) # when to start plotting lines back in time
            if istart<0: istart=0

            if i==0:

                # Plot drifter locations
                ax.plot(xp[ind20,i].T, yp[ind20,i].T, 'o', color=rgb[0,:], ms=0.6, mec='None')
                ax.plot(xp[ind50,i].T, yp[ind50,i].T, 'o', color=rgb[1,:], ms=0.6, mec='None')
                ax.plot(xp[ind100,i].T, yp[ind100,i].T, 'o', color=rgb[2,:], ms=0.6, mec='None')
                ax.plot(xp[ind500,i].T, yp[ind500,i].T, 'o', color=rgb[3,:], ms=0.6, mec='None')
                ax.plot(xp[ind3500,i].T, yp[ind3500,i].T, 'o', color=rgb[4,:], ms=0.6, mec='None')

            else:

                # Plot drifter locations
                ax.plot(xp[ind20,istart:i+1].T, yp[ind20,istart:i+1].T, '-', color=rgb[0,:], lw=0.5, zorder=0)
                ax.plot(xp[ind50,istart:i+1].T, yp[ind50,istart:i+1].T, '-', color=rgb[1,:], lw=0.5, zorder=0)
                ax.plot(xp[ind100,istart:i+1].T, yp[ind100,istart:i+1].T, '-', color=rgb[2,:], lw=0.5, zorder=0)
                ax.plot(xp[ind500,istart:i+1].T, yp[ind500,istart:i+1].T, '-', color=rgb[3,:], lw=0.5, zorder=0)
                ax.plot(xp[ind3500,istart:i+1].T, yp[ind3500,istart:i+1].T, '-', color=rgb[4,:], lw=0.5, zorder=0)


        else:

            # Plot drifter locations
            ax.plot(xp[ind20,i].T, yp[ind20,i].T, 'o', color=rgb[0,:], ms=ms, mec='None', zorder=0)
            ax.plot(xp[ind50,i].T, yp[ind50,i].T, 'o', color=rgb[1,:], ms=ms, mec='None', zorder=0)
            ax.plot(xp[ind100,i].T, yp[ind100,i].T, 'o', color=rgb[2,:], ms=ms, mec='None', zorder=0)
            ax.plot(xp[ind500,i].T, yp[ind500,i].T, 'o', color=rgb[3,:], ms=ms, mec='None', zorder=0)
            ax.plot(xp[ind3500,i].T, yp[ind3500,i].T, 'o', color=rgb[4,:], ms=ms, mec='None', zorder=0)

        if docontour:
            # Overlay surface salinity
            ax.contour(grid.x_rho.T, grid.y_rho.T, salt, [33], colors='0.1', zorder=12, linewidths=2, alpha=0.7)

    # # Time
    # ax.text(0.075, 0.95, dates[i].isoformat()[:-6], transform=ax.transAxes, fontsize=20)

    if doriver:
        # Mississippi river discharge rate
        axr = fig.add_axes([0.35, 0.2, 0.6, .2])
        # axr.set_frame_on(False) # kind of like it without the box
        for axis in ['top','left','right']:
            axr.spines[axis].set_linewidth(0.05)
        axr.spines['bottom'].set_linewidth(0.0)
        # make background rectangle so lines don't overlap
        axr.fill_between(tRiver[itstartRiver:itriver+1], R[itstartRiver:itriver+1], alpha=0.5, facecolor='0.4', edgecolor='0.4', zorder=2)
        axr.plot(tRiver[itstartRiver:itriver], R[itstartRiver:itriver], '-', color='0.4')
        axr.plot(tRiver[itriver:itendRiver+1], R[itriver:itendRiver+1], '-', color='0.4', alpha=0.3)
        axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [5, 5], '-', color='0.6', lw=0.5, alpha=0.5)
        axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [10000, 10000], '-', color='0.6', lw=0.5, alpha=0.5)
        axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [20000, 20000], '-', color='0.6', lw=0.5, alpha=0.5)
        axr.plot([tRiver[itstartRiver], tRiver[itendRiver]], [30000, 30000], '-', color='0.6', lw=0.5, alpha=0.5)
        # this makes sure alignment stays consistent in different years
        axr.autoscale(axis='x', tight=True) 
        axr.set_ylim(-1000,45000) 
        # labels
        axr.text(tRiver[mticks[-3]]+16.5, 5, '0', fontsize=9, color='0.4', alpha=0.7)
        axr.text(tRiver[mticks[-3]]+16.5, 10000, '10', fontsize=9, color='0.4', alpha=0.7)
        axr.text(tRiver[mticks[-3]]+16.5, 20000, '20', fontsize=9, color='0.4', alpha=0.7)
        axr.text(tRiver[mticks[-3]]+15, 30000, r'$30\times10^3$ m$^3$s$^{-1}$', fontsize=9, color='0.4', alpha=0.7)
        axr.text(tRiver[mticks[-7]]+15, 30000, 'Mississippi discharge', fontsize=10, color='0.2')
        # ticks
        axr.get_yaxis().set_visible(False)
        axr.get_xaxis().set_visible(False)
        # label month ticks
        for i in xrange(len(mticks)):
            axr.text(tRiver[mticks[i]], 2500, mticknames[i], fontsize=9, color='0.2')
        axr.add_patch( patches.Rectangle( (0.3, 0.162), 0.7, 0.2, transform=ax.transAxes, color='white', zorder=1))    

    if docurrents:
        # Surface currents over domain, use psi grid for common locations
        u = op.resize(np.squeeze(m.variables['u'][itmodel,-1,:,:]), 0)
        v = op.resize(np.squeeze(m.variables['v'][itmodel,-1,:,:]), 1)
        u, v = rot2d(u, v, op.resize(op.resize(anglev, 0), 1))
        Q = ax.quiver(x_psi[cdy::cdy,cdx::cdx], y_psi[cdy::cdy,cdx::cdx], u[cdy::cdy,cdx::cdx], v[cdy::cdy,cdx::cdx], 
                color='k', alpha=0.8, pivot='middle', scale=40, width=0.001, zorder=1)
        qk = ax.quiverkey(Q, 0.1, 0.795, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='0.2', fontproperties={'size': '10'})

    if dowind:
        # Wind over the domain
        Uwind = w.variables['Uwind'][itwind,:,:]
        Vwind = w.variables['Vwind'][itwind,:,:]
        Uwind, Vwind = rot2d(Uwind, Vwind, anglev)
        Q = ax.quiver(xr[wdy/2::wdy,wdx::wdx], yr[wdy/2::wdy,wdx::wdx], Uwind[wdy/2::wdy,wdx::wdx], Vwind[wdy/2::wdy,wdx::wdx], 
                color='k', alpha=0.5, scale=300, pivot='middle', headlength=3, headaxislength=2.8, zorder=1)
        qk = ax.quiverkey(Q, 0.1, 0.845, 10, r'10 m$\cdot$s$^{-1}$ wind', labelcolor='0.2', fontproperties={'size': '10'})

    if docontour:
        # Overlay 100 meter isobath
        ax.contour(xr, yr, grid.h.T, [100], colors='k', alpha=0.5, linewidths=1)

    # Drifter legend
    if dotails:
        ax.plot(0.0895, 0.9, 'or', ms=10, transform=ax.transAxes) # drifter head
        ax.plot([0.075, 0.1], [0.875, 0.875], '0.3', lw=3, transform=ax.transAxes) # drifter tail #1
        ax.plot([0.075, 0.1], [0.85, 0.85], '0.5', lw=2, transform=ax.transAxes) # drifter tail #2
        ax.text(0.125, 0.89, 'Drifter location', color='r', transform=ax.transAxes, fontsize=16)
        ax.text(0.125, 0.866, '2 days prior', color='0.3', transform=ax.transAxes, fontsize=16)
        ax.text(0.125, 0.842, '5 days prior', color='0.5', transform=ax.transAxes, fontsize=16)
    else:
        cax = fig.add_axes([0.09, 0.91, 0.35, 0.025]) #colorbar axes
        # cax = fig.add_axes([0.2, 0.81, 0.15, 0.02])
        cax.imshow(gradient, aspect='auto', interpolation='none', cmap=cmapnew)
        cax.tick_params(axis='y', labelleft=False, left=False, right=False)
        cax.tick_params(axis='x', top=False, bottom=False, labelsize=14, color='0.2', labelcolor='0.2')
        cax.set_xticks(np.arange(-0.5, 5, 1.0))
        cax.set_xticklabels(('0', '20', '50', '100', '500', '3500'))
        cax.set_xlabel('Initial drifter depth [m]', fontsize=14, color='0.2')
        # legend for contour
        ax.plot([0.05, 0.075], [0.775, 0.775], '0.1', lw=2, transform=ax.transAxes, alpha=0.7)
        ax.text(0.085, 0.767, r'33 g$\cdot$kg$^{-1}\!$ isohaline', color='0.2', transform=ax.transAxes, fontsize=10)
        # ax.plot([0.075, 0.1], [0.81, 0.81], '0.1', lw=2, transform=ax.transAxes)
        # ax.text(0.11, 0.802, '33 salinity contour', color='0.1', transform=ax.transAxes, fontsize=16)


    # Update indices
    i5daysago += 5
    i2daysago += 5

    fig.savefig(fname, bbox_inches='tight', dpi=300)

    plt.close()
