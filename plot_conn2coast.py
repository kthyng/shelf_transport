'''
Present the likelihood that oil spilled in an area will reach the coastline.
Use analysis already calculated in find_coastal_path_connectivity.
Also look at the vulnerability of the coastline for the boxes.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
import tracpy.calcs
from datetime import datetime, timedelta
import glob
import op
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar
import matplotlib.patches as Patches
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset


# of the drifters that at some point enter the coastal boxes
# When do they first enter a box? Integrate together for likelihood map
# How many drifters enter each box and at what original time? Save original location


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

# load in initial drifter starting locations in grid space
d = np.load('calcs/xyp0.npz')
xp0 = d['xp0']; yp0 = d['yp0']
d.close()

bins = (100,100)
# Calculate xrange and yrange for histograms
Xrange = [grid['xpsi'].min(), grid['xpsi'].max()]
Yrange = [grid['ypsi'].min(), grid['ypsi'].max()]

# Save a histogram of the number of drifters that started in each bin
Hstartfile = 'calcs/coastconn/likelihood/Hstart.npz'
if not os.path.exists(Hstartfile):
    Hstart, xe, ye = np.histogram2d(xp0, yp0, bins=bins, 
                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
    np.savez(Hstartfile, Hstart=Hstart, xe=xe, ye=ye)
else:
    d = np.load(Hstartfile)
    Hstart = d['Hstart']
    xe = d['xe']; ye = d['ye']
    d.close()


# # Find indices of all drifters that start in the coastal boxes
# # start indices of drifters that start in each box path
pts = np.load('calcs/alongcoastconn/inds-in-coast-paths.npz')['pts']
# pt = [] # aggregated indices of drifters that start in coast boxes
# [pt.extend(pts[j]) for j in xrange(len(pts))]
# xp0coast = xp0[pt]; yp0coast = yp0[pt]


def plot_seasonal():
    '''
    Plot seasonal comparison of likelihood, either overall or just certain parts.
    '''

    cmap = 'YlGn'
    log = False
    zoomed = True # True to add in a magnified region, for the 30 days advection timing

    d = np.load('calcs/coastpaths.npz')
    pathsxy = d['pathsxy']
    d.close()

    ## Read in files ##
    filename = 'calcs/coastconn/likelihood/hist-seasonal.npz'
    if not os.path.exists(filename):
        Files = []
        Files.append(glob.glob('calcs/coastconn/likelihood/hist-20??-0[1,2].npz'))
        Files.append(glob.glob('calcs/coastconn/likelihood/hist-20??-0[7,8].npz'))
        H = np.zeros((2,6,100,100))
        ndbox = np.zeros((2,6,342))
        for i,files in enumerate(Files): # winter and summer
            numfiles = 0
            for File in files: # months/years within winter or summer
                print File
                d = np.load(File)
                Htemp = d['Hall'] # 6x100x100, overall histogram
                # if np.isnan(H[i,:,:]).sum()>1:
                #     pdb.set_trace()
                numfiles += d['numfiles']
                H[i,:,:,:] += Htemp
                days = d['days']
                xe = d['xe']; ye = d['ye']
                # pdb.set_trace()
                ndbox[i,:,:] += d['ndbox']
                d.close()
            # pdb.set_trace()
            # Divide by number of starting drifters
            H[i,:,:,:] /= (numfiles*Hstart)
        np.savez(filename, H=H, days=days, xe=xe, ye=ye, ndbox=ndbox)
    else:
        d = np.load(filename)
        H = d['H']; days = d['days']; xe = d['xe']; ye = d['ye']; ndbox = d['ndbox']
        d.close()
    ####

    XE, YE = np.meshgrid(op.resize(xe, 0), op.resize(ye, 0))

    # Loop through advection days
    for j in xrange(H.shape[1]):

        ## Plot setup ##
        fig, axarr = plt.subplots(1,2)#, sharex=True)
        fig.set_size_inches(13, 6.6125)
        fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
        for i, ax in enumerate(axarr):
           # Titles for subplots
            if i==0:
                tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
                ax.set_title('Winter')
                # ax.set_ylabel('Source box number')
                # ax.set_xlabel('Destination box number')
                # ax.pcolormesh(mat[i,:,:], cmap=cmap)
            elif i==1:
                tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
                ax.set_title('Summer')
                # ax.xaxis.set_ticklabels([])
                # ax.yaxis.set_ticklabels([])
                # ax.get_xaxis().set_visible(False)
            levels = np.linspace(0,100,11)
            # pdb.set_trace()
            if log:
                mappable = ax.contourf(XE, YE, H[i,j,:,:].T*100., cmap=cmap, levels=levels, norm=colors.LogNorm())
            else:
                mappable = ax.contourf(XE, YE, H[i,j,:,:].T*100., cmap=cmap, levels=levels)

                # Add on vulnerability of coastline
                # Need to plot the coast boxes as patches and color them according to vulnerability level
                # http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
                # we need to normalize the data to 0..1 for the full
                # range of the colormap
                fracs = ndbox[i,j,:].astype(float)/ndbox[i,j,:].max()
                norm = colors.normalize(fracs.min(), fracs.max())
    
                # Save patches together
                patches = []
                for path in pathsxy:
                    patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))
                    # ax.add_patch(patch)

                # assign shades of colormap to the patches according to values, and plot
                for thisfrac, thispatch in zip(fracs, patches):
                    color = cm.Oranges(norm(thisfrac))
                    thispatch.set_facecolor(color)
                    ax.add_patch(thispatch)

                if zoomed and j==H.shape[1]-1: # magnification for longest advection time available

                    # pdb.set_trace()

                    # Save patches together
                    patches = []
                    for path in pathsxy:
                        patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, zorder=10, edgecolor=None))
                        # ax.add_patch(patch)

                    # Inset image
                    axins = zoomed_inset_axes(ax, 2.0, loc=4) # zoom=6
                    tracpy.plotting.background(grid=grid, ax=axins, merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])
                    axins.contourf(XE, YE, H[i,j,:,:].T*100., cmap=cmap, levels=levels)
                    # assign shades of colormap to the patches according to values, and plot
                    for thisfrac, thispatch in zip(fracs, patches):
                        color = cm.Oranges(norm(thisfrac))
                        thispatch.set_facecolor(color)
                        axins.add_patch(thispatch)

                    # pdb.set_trace()

                    # subregion of the original image
                    x1, x2, y1, y2 = 86000, 340800, 465000, 715000
                    axins.set_xlim(x1,x2)
                    axins.set_ylim(y1,y2)
                    plt.xticks(visible=False)
                    plt.yticks(visible=False)
                    plt.setp(axins,xticks=[],yticks=[])
                    # draw a bbox of the region of the inset axes in the parent axes and
                    # connecting lines between the bbox and the inset axes area
                    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
                    plt.draw()
                    plt.show()

    
                
            ax.set_frame_on(False)
        cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
        if log:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            fig.savefig('figures/coastconn/likelihood/seasonal-log-' + str(days[j]) + 'days.png', bbox_inches='tight')
        else:
            cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
            cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
            fig.savefig('figures/coastconn/likelihood/seasonal-' + str(days[j]) + 'days.png', bbox_inches='tight')
        ####
    



def likelihood():
    '''
    Aggregate likelihood of connection from locations with coast 
    in different time periods.
    '''

    nd = np.load('calcs/xyg0.npz')['xg0'].size # # of drifters

    # Loop through along-coast boxes to find which other boxes they are connected to
    years = np.arange(2013,2015)
    months = [1,2,7,8]
    days = np.array([3,5,10,15,20,30])
    for year in years:
        for month in months:
            fname = 'calcs/coastconn/likelihood/hist-' + str(year) + '-' + str(month).zfill(2) + '.npz'
            if not os.path.exists(fname):
                Files = glob.glob('calcs/alongcoastconn/' + str(year) \
                            + '-' + str(month).zfill(2) + '-*T0*.npz')

                # likelihood histogram for each advection time examined
                Hall = np.zeros((days.size, bins[0], bins[1]))
                # number of drifters reaching each coast box (yes or no, no weighting) in day days
                ndbox = np.zeros((days.size, len(pts)))
                # Histograms of where drifters originate for advection times by coast box
                H = np.zeros((days.size, len(pts), bins[0], bins[1]))
                
                for File in Files:
                    print File
                    d = np.load(File)
                    # [# of box paths x # drifters that enter a box x 5 (max # of crosses checked for)]
                    inbox = d['inbox'] # time in decimal days when a drifter enters a box path
                    # outbox = d['outbox'] # time in decimal days when a drifter exists a box path
                    # inds are referenced to the drifters in the shelf transport runs
                    inds = d['iinside'] # indices from the original drifters corresponding to in/outbox
                    d.close()

                    # code to switch between sets of indices
                    # this has Trues for the drifters for this simulation that enter 
                    # the outer path
                    code = np.zeros(nd); code[inds] = 1; code = find(code.astype(bool))
                    # xp0[code] gives the x positions in projected space of drifters that were
                    # analyzed in the code due to entering the coastal boxes at some point.

                    # xp0temp = xp0[code]
                    # xp0temp = xp0temp[np.newaxis,:].repeat(inbox.shape[0], axis=0)
                    # yp0temp = yp0[code]
                    # yp0temp = yp0temp[np.newaxis,:].repeat(inbox.shape[0], axis=0)

                    for i, day in enumerate(days): # loop through number of advection days

                        # Drifters that enter a coast box within day days [coast box x set of drifters]
                        ind = (inbox[:,:,0]<=day)

                        # How many drifters enter each box by day days?
                        ndbox[i,:] += ind.sum(axis=1)

                        # overall histogram
                        # which drifters enter a coast box in this time period
                        tocoast = (ind.sum(axis=0)).astype(bool)
                        # xp = xp0[code][ind]; yp = yp0[code][ind]
                        # # Original xp, yp locations of drifters that enter a coast box within day days
                        # xpsaveall = xp0temp[ind]; ypsaveall = yp0temp[ind]
                        xpsaveall = xp0[code][tocoast]; ypsaveall = yp0[code][tocoast];
                        Halltemp, _, _ = np.histogram2d(xpsaveall, ypsaveall, bins=bins, 
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                        Hall[i,:,:] += Halltemp
                        # pdb.set_trace()

                        # loop through each coast box path to calculate separate origin histogram
                        for j in xrange(len(pts)):

                            # indices of drifters that start in this box, referenced to shelf transport seeding
                            pt = pts[j] 

                            # projected drifter origin locations that end up in this coast box in day days
                            xpsave = xp0[code][ind[j,:]]; ypsave = yp0[code][ind[j,:]]

                            # pdb.set_trace()

                            # Add in drifters that start in coast boxes
                            xpcoast = set(xp0[pt]) - set(xpsave)
                            xpsave = np.concatenate((xpsave, list(xpcoast)))
                            ypcoast = set(yp0[pt]) - set(ypsave)
                            ypsave = np.concatenate((ypsave, list(ypcoast)))

                            # Find histogram of xpsave, ypsave points for this simulation/box origin points
                            Htemp, _, _ = np.histogram2d(xpsave, ypsave, bins=bins, 
                                            range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                            # aggregate number of drifters starting in different histogram bins 
                            # that reach coastline for each month/year combination
                            H[i,j,:,:] += Htemp 



                # Save the month/year's worth of histograms
                # numfiles is to calculate the number drifters from bins for the the number of runs
                # aggregated together, compared with the appropriate number of starting drifters overall
                np.savez(fname, H=H, xe=xe, ye=ye, days=days, numfiles=len(Files), ndbox=ndbox, Hall=Hall)
                # pdb.set_trace()
            
        

if __name__ == "__main__":
    likelihood()     

