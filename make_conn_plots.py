'''
Read in connectivity calculations from find_coastal_path_connectivity.py
and make plots.
'''

import matplotlib as mpl
# mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
# import tracpy
# import tracpy.plotting
# import tracpy.calcs
from datetime import datetime, timedelta
import glob
# import op
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar
import cmocean.cm as cmo
import matplotlib.patches as Patches
import pandas as pd


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
colu = '#218983'  # upcoast color
cold = '#cb6863'  # downcoast color


# # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# loc = '../grid.nc'
# proj = tracpy.tools.make_proj('nwgom')
# grid = tracpy.inout.readgrid(loc, proj)

# load in paths for coastline boxes
d = np.load('calcs/coastpaths.npz', encoding='latin1') # use paths in grid space
pathsxy = d['pathsxy']
outerpathxy = d['outerpathxy'].item()
d.close()

# distance along the coast boxes
# code from plot_sampledrifters.py
dist = np.zeros(len(pathsxy))
verts0 = pathsxy[0].vertices
for i, path in enumerate(pathsxy):
    verts1 = path.vertices
    dist[i:] += np.sqrt((verts1[0,0]-verts0[0,0])**2+(verts1[0,1]-verts0[0,1])**2)
    verts0 = verts1.copy()
dist /= 1000 # convert to km
dmax = dist.max()

X, Y = np.meshgrid(dist, dist)

# dic = {}
# Port Mansfield, Port Aransas, Port O'Connor, Galveston, Atchafalaya,
# Terrebonne, Barataria
names = ['bpm', 'bpa', 'bpoc', 'bgalv', 'batch', 'bterr', 'bbara', 'all']
distances = [400, 555, 645, 840, 1175, 1280, 1350]  # atch 1190, bterr1290
boxes = [np.arange(76,88), np.arange(107,119), np.arange(126,138), np.arange(165,177),
         np.arange(234,246), np.arange(257,269), np.arange(271, 283), np.arange(0,342)]
# # for name, dist, box in zip(names, dists, boxes):
# for name, distance in zip(names, distances):
#     dic[name] = {}
#     dic[name]['dist'] = distance
#     # dic[name]['boxes'] = box

boxdict = {'bpm': np.arange(76,88), 'bpa': np.arange(107,119),
         'bpoc': np.arange(126,138), 'bgalv': np.arange(165,177),
         'batch': np.arange(234,246), 'bterr': np.arange(257,269),
         'bbara': np.arange(271, 283), 'all': np.arange(0,342) }


def plot_domain():
    '''
    Plot explanatory map for along-shore box key.
    '''

    fig = plt.figure()
    ax = fig.add_subplot(111)
    tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))

    # Add outer path
    ax.plot(outerpathxy.vertices[:,0], outerpathxy.vertices[:,1], 'b')

    # Plot labels for boxes
    # if plottype == 'dist':
    for (i, path), dis in zip(enumerate(pathsxy), dist):
        # if np.mod(i, 10) == 0:
        #     ax.text(path.vertices[0,0], path.vertices[0,1], str(dis), fontsize=12, color='k')
            # ax.text(path.vertices[0,0], path.vertices[0,1], str(i), fontsize=15, color='g')
        if (abs(dis-np.asarray(distances)) < 2.5).any():
            ind = np.where(abs(dis-np.asarray(distances)) < 2.5)[0]
            ax.text(path.vertices[0,0], path.vertices[0,1], str(distances[ind]), fontsize=20, color='r')
            # ax.text(path.vertices[0,0], path.vertices[0,1], str(i), fontsize=15, color='g')

            patches = []
            for j in boxes[ind]:
                ax.add_patch(Patches.PathPatch(pathsxy[j], facecolor='yellow'))
                # if (j == boxes[ind][0]) or (j == boxes[ind][-1]):
                #     ax.text(path.vertices[0,0], path.vertices[0,1], str(j), fontsize=15, color='g')
            #     patches.append(Patches.PathPatch(pathxsy[j]))#, facecolor='orange', lw=0, zorder=10, edgecolor=None))
            # # assign shades of colormap to the patches according to values, and plot
            # for thisfrac, thispatch in zip(fracs, patches):
            #     color = cmo.matter(norm(thisfrac))
            #     thispatch.set_facecolor(color)
            #     ax.add_patch(thispatch)

    # elif plottype == 'boxnum':
    # for i, path in enumerate(pathsxy):
    #     if np.mod(i, 10) == 0:
    #         ax.text(path.vertices[0,0], path.vertices[0,1], str(i), fontsize=15, color='r')
        # # only plot every 50th
        # if np.mod(i,50)==0:
        #     if i<110:
        #         dx = 11000; dy = 0; # shift
        #     elif i==150:
        #         dx = 5000; dy = -30000
        #     elif i>150 and i<300:
        #         dx = 0; dy = -50000
        #     elif i==300:
        #         dy = -40000
        #     ax.text(path.vertices[0,0] + dx, path.vertices[0,1] + dy, str(i), fontsize=15, color='r')

    # fig.savefig('figures/alongcoastconn/domain.png', bbox_inches='tight')


def plot_bays_seasonal():
    '''Plot starting and ending for bays'''

    cmap = cmo.curl_r

    filename = 'calcs/alongcoastconn/conn-seasonal.npz'
    mat = np.load(filename)['mat']  # 2x342x342
    mat = mat.transpose([0,2,1])  # transpose needed since pcolormesh flips it?
    bmax = abs(mat[0,boxes[10],boxes[10]].sum())
    dbox = 14

    fig = plt.figure(figsize=(8.5,11))
    fig.subplots_adjust(wspace=0.05)
    # winter
    ax = fig.add_subplot(1,2,1)
    ax.set_title('Winter')
    i = 0
    for box in np.arange(dbox, len(dist)-dbox, dbox):

        # traveling to: downcoast
        y = (mat[i,:box,box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[:box], y+dist[box], dist[box], color=cold)
        # traveling to: upcoast
        y = (mat[i,box+1:,box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[box+1:], y+dist[box], dist[box], color=colu)

        # use one box earlier
        box -= 1
        # traveling from: upcoast
        y = -(mat[i,box,:box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[:box], y+dist[box], dist[box], color=colu)
        # traveling from: downcoast
        y = -(mat[i,box,box+1:]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[box+1:], y+dist[box], dist[box], color=cold)
        # # zero
        # ax.plot(dist, np.zeros(dist.size)+dist[box], lw=0.5, color='0.1', alpha=0.5)
    ax.axis('tight')
    ax.set_xticks(np.arange(0, dist.max(), 500))
    ax.set_xticks(np.arange(0, dist.max(), 100), minor=True)
    ax.set_yticks(np.arange(0, dist.max(), 100), minor=True)
    ax.set_yticklabels([''])
    ax.grid(lw=0.1, color='0.2', linestyle='-', alpha=0.6, which='both')
    ax.set_frame_on(False)
    ax.set_xlabel('along-coast distance [km]')
    # summer
    ax = fig.add_subplot(1,2,2)
    ax.set_title('Summer')
    i = 1
    for box in np.arange(dbox, len(dist)-dbox, dbox):

        # traveling to: downcoast
        y = (mat[i,:box,box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[:box], y+dist[box], dist[box], color=cold)
        # traveling to: upcoast
        y = (mat[i,box+1:,box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[box+1:], y+dist[box], dist[box], color=colu)

        # use one box earlier
        box -= 1
        # traveling from: upcoast
        y = -(mat[i,box,:box]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[:box], y+dist[box], dist[box], color=colu)
        # traveling from: downcoast
        y = -(mat[i,box,box+1:]/bmax)*40
        y[y==0] = np.nan  # don't want to include 0 values
        ax.fill_between(dist[box+1:], y+dist[box], dist[box], color=cold)
        # # zero
        # ax.plot(dist, np.zeros(dist.size)+dist[box], lw=0.5, color='0.1', alpha=0.5)
    ax.axis('tight')
    ax.set_xticks(np.arange(0, dist.max(), 500))
    ax.set_xticks(np.arange(0, dist.max(), 100), minor=True)
    ax.set_yticks(np.arange(0, dist.max(), 100), minor=True)
    ax.set_yticklabels([''])
    ax.grid(lw=0.1, color='0.2', linestyle='-', alpha=0.6, which='both')
    ax.set_frame_on(False)
    ax.set_xlabel('along-coast distance [km]')
    fig.savefig('figures/alongcoastconn/bays.png', bbox_inches='tight')



def plot_interannual():
    '''
    Use calculated files from run() to plot connectivity matrix.
    '''

    cmap = cmo.curl_r
    season = 'winter'  # 'winter' or 'summer'
    log = False
    regions = True  # do plot with region lines (MX, TX, LA)
    baylabels = False  # mark locations with ticks
    bayvalues = False  # annotate bay values
    largefonts = False  # use large fonts for presentation plot

    xticklocs = np.arange(0, 2000, 500)
    years = np.arange(2004, 2015)

    if season == 'winter':
        iseason = 0
    elif season == 'summer':
        iseason = 1

    ## Read in files ##
    filename = 'calcs/alongcoastconn/conn-interannual.npz'
    if not os.path.exists(filename):
        mat = np.zeros((years.size, 2,342,342))
        for j, year in enumerate(years):
            Files = []
            Files.append(glob.glob('calcs/alongcoastconn/conn-' + str(year) + '-0[1,2].npz'))
            Files.append(glob.glob('calcs/alongcoastconn/conn-' + str(year) + '-0[7,8].npz'))
            for i,files in enumerate(Files): # winter and summer
                for File in files: # months/years within winter or summer
                    print(File)
                    d = np.load(File)
                    mat[j, i,:,:] += d['mat']
                    if np.isnan(mat[j, i,:,:]).sum()>1:
                        pdb.set_trace()
                mat[j, i,:,:] /= len(files)
        np.savez(filename, mat=mat)
    else:
        mat = np.load(filename)['mat']
    ####

    # make one side of triangle positive and other negative
    ix, iy = np.tril_indices(mat.shape[2], k=1)
    # mat[0,:,:] = -np.triu(mat[0,:,:].T, k=1)
    # mat[1,:,:] = -np.triu(mat[1,:,:].T, k=1)
    mat[:, :, ix, iy] = -mat[:, :, ix, iy]

    fig, axarr = plt.subplots(4,3)
    fig.set_size_inches(8.9, 13.5)
    fig.subplots_adjust(left=0.04, bottom=0.1, right=1.0, top=0.99, wspace=0.08, hspace=0.1)

    for i, ax in enumerate(axarr.flat):
        ax.set_frame_on(False)
        if i == 11:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            continue
        ax.set_title(str(years[i]))
        if i == 9:
            ax.set_ylabel('Along coast start location [km]')
            ax.set_xlabel('Along coast end location [km]')
        ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
        ax.yaxis.set_ticks_position('left')
        ax.axis('equal')
        ax.set_xticks(xticklocs);
        ax.set_yticks(xticklocs);
        if i != 9:
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])
        if log:
            mappable = ax.pcolormesh(mat[i,iseason,:,:]*100., cmap=cmap, vmin=1., vmax=100., norm=colors.LogNorm())
        else:
            mappable = ax.pcolormesh(X, Y, mat[i,iseason,:,:]*100., cmap=cmap, vmax=100., vmin=-100)
        # import pdb; pdb.set_trace()
        # Plot a few ticks for notable locations
        left = -25; right = 35
        for distance, box in zip(distances, boxes):
            # max value for normalization is the boxes to themselves
            bmax = abs(mat[i,iseason,box,box].sum())
            if baylabels:
                ax.autoscale(enable=False)
                # horizontal
                ax.plot([left, right], [distance, distance], '-', color='0.3', lw=0.5, alpha=0.7)
                # vertical
                ax.plot([distance, distance], [left, right], '-', color='0.3', lw=0.5, alpha=0.7)
            if bayvalues:
                # plot Connectivity for all other bays
                for distance2, box2 in zip(distances, boxes):
                    if distance2 == distance:
                        continue
                    boxval = (mat[i,iseason,box,box2].sum()/bmax)*100
                    if not boxval == 0:
                        # pass
                        # ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                        #            cmap=cmo.curl_r, vmin=-100, vmax=100,
                        #            linewidths=0.1)
                        # print('not zero:', boxval)
                        if largefonts:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=17,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                        else:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=9,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                    else:
                        # Otherwise the axes get moved by scatter
                        # http://stackoverflow.com/questions/19916295/pyplot-scatter-changes-the-data-limits-of-the-axis
                        ax.autoscale(enable=False)
                        if largefonts:
                            ax.scatter(distance2, distance, s=70, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        else:
                            ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        # print('zero:', boxval)

        if regions:
            # Overlay lines boxes for region of coastline
            # horizontal: Mexico-Texas border
            ax.plot([0, dmax], [340, 340], '-', color='k', lw=2, alpha=0.1)
            # horizontal: Texas-Louisiana border
            ax.plot([0, dmax], [940, 940], '-', color='k', lw=2, alpha=0.1)
            # vertical: Mexico-Texas border
            ax.plot([340, 340], [0, dmax], '-', color='k', lw=2, alpha=0.1)
            # vertical: Texas-Louisiana border
            ax.plot([940, 940], [0, dmax], '-', color='k', lw=2, alpha=0.1)
    # cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
    cax = fig.add_axes([0.25, 0.03, 0.5, 0.02]) #colorbar axes
    ticklabels = ['100', '80', '60', '40', '20', '0', '20', '40', '60', '80', '100']
    if log:
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
        cb.set_label('Connectivity [%]')
        fig.savefig('figures/alongcoastconn/interannual-log.png', bbox_inches='tight')
    else:
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')#, pad=0.18)
        cb.set_label('Connectivity [%]')
        cb.set_ticks(np.arange(-100, 120, 20))
        cb.set_ticklabels(ticklabels)
        fname = 'figures/alongcoastconn/interannual-' + season
        if regions:
            fname += '-regions'
        if baylabels:
            fname += '-baylabels'
        if bayvalues:
            fname += '-bayvalues'
        if largefonts:
            fname += '-largefonts'
        fig.savefig(fname + '.png', bbox_inches='tight', dpi=300)
    ####


def plot_seasonal():
    '''
    Use calculated files from run() to plot connectivity matrix.
    '''

    cmap = cmo.curl_r  # 'YlGn'
    log = False
    regions = True  # do plot with region lines (MX, TX, LA)
    baylabels = True  # mark locations with ticks
    bayvalues = True  # annotate bay values
    largefonts = True  # use large fonts for presentation plot

    xticklocs = np.arange(0, 2000, 500)

    ## Read in files ##
    filename = 'calcs/alongcoastconn/conn-seasonal.npz'
    if not os.path.exists(filename):
        Files = []
        Files.append(glob.glob('calcs/alongcoastconn/conn-20??-0[1,2].npz'))
        Files.append(glob.glob('calcs/alongcoastconn/conn-20??-0[7,8].npz'))
        mat = np.zeros((2,342,342))
        for i,files in enumerate(Files): # winter and summer
            for File in files: # months/years within winter or summer
                print(File)
                d = np.load(File)
                mat[i,:,:] += d['mat']
                if np.isnan(mat[i,:,:]).sum()>1:
                    pdb.set_trace()
            mat[i,:,:] /= len(files)
        np.savez(filename, mat=mat)
    else:
        mat = np.load(filename)['mat']
    ####

    # make one side of triangle positive and other negative
    ix, iy = np.tril_indices(mat.shape[1], k=1)
    # mat[0,:,:] = -np.triu(mat[0,:,:].T, k=1)
    # mat[1,:,:] = -np.triu(mat[1,:,:].T, k=1)
    mat[:, ix, iy] = -mat[:, ix, iy]

    ## Plot setup ##
    fig, axarr = plt.subplots(1,2)#, sharex=True)
    fig.set_size_inches(11.0, 6.6125)
    fig.subplots_adjust(left=0.045, bottom=0.175, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
    for i, ax in enumerate(axarr):
       # Titles for subplots
        if i==0:
            # tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title('Winter')
            ax.set_ylabel('Along coast start location [km]')
            ax.set_xlabel('Along coast end location [km]')
            ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
            ax.yaxis.set_ticks_position('left')
            ax.axis('equal')
            ax.set_xticks(xticklocs);
            ax.set_yticks(xticklocs);
            # ax.pcolormesh(mat[i,:,:], cmap=cmap)
        elif i==1:
            # tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            ax.set_xlabel('Along coast end location [km]')
            ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
            ax.yaxis.set_ticks_position('left')
            ax.axis('equal')
            # ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.set_xticks(xticklocs);
            ax.set_yticks(xticklocs);
            # ax.get_xaxis().set_visible(False)
        if log:
            mappable = ax.pcolormesh(mat[i,:,:]*100., cmap=cmap, vmin=1., vmax=100., norm=colors.LogNorm())
        else:
            mappable = ax.pcolormesh(X, Y, mat[i,:,:]*100., cmap=cmap, vmax=100., vmin=-100)
        # import pdb; pdb.set_trace()
        ax.set_frame_on(False)
        # Plot a few ticks for notable locations
        left = -25; right = 35
        for distance, box in zip(distances, boxes):
            # max value for normalization is the boxes to themselves
            bmax = abs(mat[i,box,box].sum())
            if baylabels:
                ax.autoscale(enable=False)
                # horizontal
                ax.plot([left, right], [distance, distance], '-', color='0.3', lw=0.5, alpha=0.7)
                # vertical
                ax.plot([distance, distance], [left, right], '-', color='0.3', lw=0.5, alpha=0.7)
            if bayvalues:
                # plot Connectivity for all other bays
                for distance2, box2 in zip(distances, boxes):
                    if distance2 == distance:
                        continue
                    boxval = (mat[i,box,box2].sum()/bmax)*100
                    if not boxval == 0:
                        # pass
                        # ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                        #            cmap=cmo.curl_r, vmin=-100, vmax=100,
                        #            linewidths=0.1)
                        # print('not zero:', boxval)
                        if largefonts:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=17,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                        else:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=9,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                    else:
                        # Otherwise the axes get moved by scatter
                        # http://stackoverflow.com/questions/19916295/pyplot-scatter-changes-the-data-limits-of-the-axis
                        ax.autoscale(enable=False)
                        if largefonts:
                            ax.scatter(distance2, distance, s=70, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        else:
                            ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        # print('zero:', boxval)

        if regions:
            # Overlay lines boxes for region of coastline
            # horizontal: Mexico-Texas border
            ax.plot([0, dmax], [340, 340], '-', color='k', lw=2, alpha=0.1)
            # horizontal: Texas-Louisiana border
            ax.plot([0, dmax], [940, 940], '-', color='k', lw=2, alpha=0.1)
            # vertical: Mexico-Texas border
            ax.plot([340, 340], [0, dmax], '-', color='k', lw=2, alpha=0.1)
            # vertical: Texas-Louisiana border
            ax.plot([940, 940], [0, dmax], '-', color='k', lw=2, alpha=0.1)
    cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
    ticklabels = ['100', '80', '60', '40', '20', '0', '20', '40', '60', '80', '100']
    # if log:
    #     cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
    #     cb.set_label('Connectivity [%]')
    #     fig.savefig('figures/alongcoastconn/seasonal-log.png', bbox_inches='tight')
    # else:
    #     cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')#, pad=0.18)
    #     cb.set_label('Connectivity [%]')
    #     cb.set_ticks(np.arange(-100, 120, 20))
    #     cb.set_ticklabels(ticklabels)
    #     fname = 'figures/alongcoastconn/seasonal'
    #     if regions:
    #         fname += '-regions'
    #     if baylabels:
    #         fname += '-baylabels'
    #     if bayvalues:
    #         fname += '-bayvalues'
    #     if largefonts:
    #         fname += '-largefonts'
    #     fig.savefig(fname + '.png', bbox_inches='tight', dpi=300)
    ####


def plot_monthly():
    '''
    Use calculated files from run() to plot connectivity matrix for each month
    in the year together.
    '''

    cmap = cmo.curl_r
    log = False
    regions = True  # do plot with region lines (MX, TX, LA)
    baylabels = False  # mark bay locations with ticks
    baylines = True  # mark bay locations with lines
    bayvalues = False  # annotate bay values
    largefonts = False  # use large fonts for presentation plot

    xticklocs = np.arange(0, 2000, 500)
    months = np.arange(1, 13)

    ## Read in files ##
    # filename = 'calcs/alongcoastconn/conn-interannual.npz'
    # if not os.path.exists(filename):
    mat = np.zeros((months.size,342,342))
    Files = []
    for j, month in enumerate(months):
        Files.append(glob.glob('calcs/alongcoastconn/conn-20??-' + str(month).zfill(2) + '.npz'))
    for j,files in enumerate(Files): # months
        for File in files: # years within month
            # print(File)
            d = np.load(File)
            mat[j, :,:] += d['mat']
            # if np.isnan(mat[j, i,:,:]).sum()>1:
            #     pdb.set_trace()
        mat[j, :,:] /= len(files)
    # np.savez(filename, mat=mat)
    # else:
    #     mat = np.load(filename)['mat']
    ####

    # make one side of triangle positive and other negative
    ix, iy = np.tril_indices(mat.shape[2], k=1)
    mat[:, ix, iy] = -mat[:, ix, iy]

    fig, axarr = plt.subplots(4,3, sharex=True, sharey=True, figsize=(8.9, 13.5))
    fig.subplots_adjust(left=0.1, bottom=0.12, right=0.98, top=0.97, wspace=0.06, hspace=0.1)

    for i, ax in enumerate(axarr.flat):
        ax.set_frame_on(False)
        ax.set_title(datetime(1970,months[i],1,0,0).strftime('%b'))
        if i == 9:
            ax.set_ylabel('Along coast start location [km]')
            ax.set_xlabel('Along coast end location [km]')
            # import pdb; pdb.set_trace()
        if log:
            mappable = ax.pcolormesh(mat[i,:,:]*100., cmap=cmap, vmin=1., vmax=100., norm=colors.LogNorm())
        else:
            mappable = ax.pcolormesh(X, Y, mat[i,:,:]*100., cmap=cmap, vmax=100., vmin=-100)
        ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
        ax.yaxis.set_ticks_position('left')
        ax.axis('equal')
        ax.set_xticks(xticklocs);
        ax.set_yticks(xticklocs);
        # if i != 9:
        #     ax.yaxis.set_ticklabels([])
        #     ax.xaxis.set_ticklabels([])
        # import pdb; pdb.set_trace()
        # Plot a few ticks for notable locations
        left = -25; right = 35
        for distance, box in zip(distances, boxes):
            # max value for normalization is the boxes to themselves
            bmax = abs(mat[i,box,box].sum())
            if baylabels:
                ax.autoscale(enable=False)
                # horizontal
                ax.plot([left, right], [distance, distance], '-', color='0.3', lw=0.5, alpha=0.7)
                # vertical
                ax.plot([distance, distance], [left, right], '-', color='0.3', lw=0.5, alpha=0.7)
            if baylines:
                ax.autoscale(enable=False)
                # horizontal
                ax.plot([0, dmax], [distance, distance], '-', color='0.3', lw=0.1, alpha=0.5)
                # vertical
                ax.plot([distance, distance], [0, dmax], '-', color='0.3', lw=0.1, alpha=0.5)
            if bayvalues:
                # plot Connectivity for all other bays
                for distance2, box2 in zip(distances, boxes):
                    if distance2 == distance:
                        continue
                    boxval = (mat[i,box,box2].sum()/bmax)*100
                    if not boxval == 0:
                        # pass
                        # ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                        #            cmap=cmo.curl_r, vmin=-100, vmax=100,
                        #            linewidths=0.1)
                        # print('not zero:', boxval)
                        if largefonts:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=17,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                        else:
                            ax.text(distance2+1, distance-16, '%d' % abs(boxval), fontsize=9,
                                    horizontalalignment='center', alpha=0.8, color='0.2')
                    else:
                        # Otherwise the axes get moved by scatter
                        # http://stackoverflow.com/questions/19916295/pyplot-scatter-changes-the-data-limits-of-the-axis
                        ax.autoscale(enable=False)
                        if largefonts:
                            ax.scatter(distance2, distance, s=70, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        else:
                            ax.scatter(distance2, distance, s=50, c=boxval, marker='s',
                                       cmap=cmo.curl_r, vmin=-100, vmax=100,
                                       linewidths=0.2)
                        # print('zero:', boxval)

        if regions:
            # Overlay lines boxes for region of coastline
            # horizontal: Mexico-Texas border
            ax.plot([0, dmax], [340, 340], '-', color='k', lw=2, alpha=0.1)
            # horizontal: Texas-Louisiana border
            ax.plot([0, dmax], [940, 940], '-', color='k', lw=2, alpha=0.1)
            # vertical: Mexico-Texas border
            ax.plot([340, 340], [0, dmax], '-', color='k', lw=2, alpha=0.1)
            # vertical: Texas-Louisiana border
            ax.plot([940, 940], [0, dmax], '-', color='k', lw=2, alpha=0.1)
    cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
    # cax = fig.add_axes([0.25, 0.03, 0.5, 0.02]) #colorbar axes
    ticklabels = ['100', '80', '60', '40', '20', '0', '20', '40', '60', '80', '100']
    if log:
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
        cb.set_label('Connectivity [%]')
        fig.savefig('figures/alongcoastconn/monthly-log.png', bbox_inches='tight')
    else:
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')#, pad=0.18)
        cb.set_label('Connectivity [%]')
        cb.set_ticks(np.arange(-100, 120, 20))
        cb.set_ticklabels(ticklabels)
        fname = 'figures/alongcoastconn/monthly'
        if regions:
            fname += '-regions'
        if baylabels:
            fname += '-baylabels'
        if baylines:
            fname += '-baylines'
        if bayvalues:
            fname += '-bayvalues'
        if largefonts:
            fname += '-largefonts'
        fig.savefig(fname + '.png', dpi=300)#, bbox_inches='tight')
    ####


def plot_bayconn(boxnameto, boxnamefrom):
    '''Plot connectivity between specific bays throughout the year.'''

    base = 'calcs/alongcoastconn/conn_in_time/'

    # indices of boxes to which drifters are traveling
    iboxto = boxdict[boxnameto]
    # indices of boxes from which drifters are traveling
    iboxfrom = boxdict[boxnamefrom]

    # load in dataframe from being calculated in calc_bayconn()
    df = pd.read_csv(base + 'to_' + boxnameto + '.csv', parse_dates=True, index_col=0)

    # calculate day of year with higher res for x axis
    df['doy'] = df.index.dayofyear + df.index.hour/24.

    import seaborn as sns
    # fig, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(7, 3))
    # fig.subplots_adjust(left=0.04, bottom=0.1, right=1.0, top=0.97, wspace=0.06, hspace=0.1)

    # for i, ax in enumerate(axarr.flat):
    # ax.set_frame_on(False)
    # (df['bpoc-30']['2004']/df['nsims']['2004']).plot()

    # calculate quintile to quintile with groupby and show with fill between
    # (df['bpoc-30']/df['nsims']).groupby(df['doy']).mean().plot()
    # calculate quantile, which has 2 results for each index for the two quantiles
    quantlow = df['bpoc-30'].groupby(df['doy']).quantile(0.367)
    quanthigh = df['bpoc-30'].groupby(df['doy']).quantile(0.733)
    # quantlow = (df['bpoc-30']/df['nsims']).groupby(df['doy']).quantile(0.367)
    # quanthigh = (df['bpoc-30']/df['nsims']).groupby(df['doy']).quantile(0.733)
    plt.fill_between(quantlow.index, quantlow, quanthigh, alpha=0.2)

    df2 = (df['bpoc-30']/df['nsims']).groupby(df['doy']).mean()
    (df['bpoc-30']/df['nsims']).groupby(df['doy']).mean().plot()
    (df['bpoc-20']/df['nsims']).groupby(df['doy']).mean().plot()

    # label ticks at start of month
    # ticks for months on river discharge
    mticks = [bisect.bisect_left(datesRiver, monthdate) for monthdate in np.asarray(monthdates)]
    mticknames = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    # how to easily get 1st day of each month in terms of day of year


def calc_bayconn():
    '''Combine together files from run_with_times() to use in plot_bayconn.

    Files are of the types calcs/alongcoastconn/conn_in_time/20??-??-??T0?.npz
    '''
    base = 'calcs/alongcoastconn/conn_in_time/'

    ndays = np.array([5, 10, 15, 20, 25, 30])  # number of advection days to consider

    # create a dataframe for transport to each area, for all of the time examined, 4 hourly
    dfdates = pd.date_range(start='2004-01-01 00:00', end='2014-10-01 00:00', freq='14400S')
    # make a dictionary of dataframes. They are sorted by where transport is going to.
    dfs = {}
    for tokey in boxdict.keys():  # dataframes are by transport "to"
        dfs[tokey] = pd.DataFrame(index=dfdates)
        for fromkey in boxdict.keys():  # each has a column of transport "from"
            for nday in ndays:
                # add in column for ndrifters
                dfs[tokey][fromkey + '-' + str(nday)] = 0
        # add in column for # simulations integrated in ndrifters column
        dfs[tokey]['nsims'] = 0

    Files = glob.glob(base + '20??-??-??T0?.npz')

    for File in Files:  # loop through all daily simulation files
        d = np.load(File)
        # mat: time x from coast boxes x to coast boxes, t: 4 hourly times
        mat = d['mat']; t = d['t']
        d.close()

        startdate = File.split('/')[-1].split('.')[0]  # start date of simulation in string form
        startdatedt = datetime.strptime(startdate, '%Y-%m-%dT%H')  # date in datetime format

        for tokey in boxdict.keys():  # loop through TO areas
            for fromkey in boxdict.keys():  # loop through FROM areas
                # import pdb; pdb.set_trace()
                isfrom = boxdict[fromkey][0]  # starting index for "from"
                iefrom = boxdict[fromkey][-1]  # ending index for "from"
                isto = boxdict[tokey][0]  # starting index for "to"
                ieto = boxdict[tokey][-1]  # end index for "to"

                for nday in ndays:
                    enddatedt = startdatedt + timedelta(days=int(nday))
                    enddate = enddatedt.strftime('%Y-%m-%dT%H')  # end date in string format
                    # dates that ndrifters in simulation cover
                    dates = pd.date_range(start=startdate, end=enddate, freq='14400S')
                    # remove first date since nothing happens then
                    dates = dates[1:]

                    # sum number of drifters across times for "from" and "to" regions
                    # then average across the number of "from" boxes
                    # ndrifters ends up being in time
                    # use up to the number of times in dates
                    ndrifters = mat[:len(dates), isfrom:iefrom, isto:ieto].sum(axis=2).sum(axis=1)/boxdict[fromkey].size
                    # add to dataframe at the relevant times
                    dfs[tokey].loc[dates, fromkey + '-' + str(nday)] = ndrifters
                # keep track of number of simulations being added together
                dfs[tokey].loc[dates, 'nsims'] += 1
            # import pdb; pdb.set_trace()

    for tokey in boxdict.keys():  # dataframes are by transport "to"
        dfs[tokey].to_csv(base + 'to_' + tokey + '.csv')


def run():
    '''
    Combine previously calculated files into month/year files.
    '''

    # load in paths for coastline boxes
    d = np.load('calcs/coastpaths.npz') # use paths in grid space
    paths = d['pathsg']
    d.close()

    # nd = np.load('calcs/xyg0.npz')['xg0'].size # # of drifters

    # start indices of drifters that start in each box path
    filename = 'calcs/alongcoastconn/inds-in-coast-paths.npz'
    # indices in pts are referenced to the drifters that start
    # in the shelf transport simulations
    if os.path.exists(filename):
        pts = np.load(filename)['pts']
    else:
        # do this once and save
        pts = []
        for i,path in enumerate(paths):
            # which points are inside the regions
            pts.append(find(path.contains_points(np.vstack((xg0.flat, yg0.flat)).T).reshape(xg0.shape)))
        np.savez(filename, pts=pts)

    # Loop through along-coast boxes to find which other boxes they are connected to
    # and allow for the 5 crossings kept track of
    mat = np.zeros((len(paths),len(paths)))  # original along-coast conn, for 30 days
    # inmat = np.zeros((len(paths),len(paths), 5))  # drifters going into boxes
    # outmat = np.zeros((len(paths),len(paths), 5))  # drifters going out of boxes
    years = np.arange(2004,2015)
    months = np.arange(1, 13)
    # months = [1,2,7,8]
    for year in years:
        for month in months:
            matfile = 'calcs/alongcoastconn/conn-' + str(year) + '-' + str(month).zfill(2) + '.npz'
            if not os.path.exists(matfile):
                Files = glob.glob('calcs/alongcoastconn/' + str(year) \
                            + '-' + str(month).zfill(2) + '-*T0*.npz')

                for File in Files:
                    print(File)
                    d = np.load(File)
                    # [# of box paths x # drifters that enter a box x 5 (max # of crosses checked for)]
                    inbox = d['inbox'] # time in decimal days when a drifter enters a box path
                    # outbox = d['outbox'] # time in decimal days when a drifter exists a box path
                    # inds are referenced to the drifters in the shelf transport runs
                    inds = d['iinside'] # indices from the original drifters corresponding to in/outbox
                    d.close()
                    # # code to switch between sets of indices
                    # # this has Trues for the drifters for this simulation that enter
                    # # the outer path
                    # code = np.zeros(nd); code[inds] = 1; code = code.astype(bool)
                    mattemp = np.eye(inbox.shape[0])
                    # For each box, how many drifters go out from it in to another box?
                    for i, path in enumerate(paths):

                        pt = pts[i] # indices of drifters that start in this box

                        # find indices of where indices in pt can be found in inds/inbox/outbox
                        # pt is the same as inds[code]
                        code = []
                        for p in pt:
                            code.extend(find(p==inds))

                        # # which paths are connected to the box in question
                        # connected = ~np.isnan(np.nansum(inbox[:,code,0], axis=1))

                        # the number of drifters that enter each box
                        # import pdb; pdb.set_trace()
                        connected = np.sum(~np.isnan(inbox[:,code,0]), axis=1)
                        # manually do the box the drifters started in
                        connected[i] = len(pt)

                        # # put a 1 for the boxes where a drifter will connect
                        # mattemp[i,connected] = 1

                        # Normalize down to 1 for drifters starting in each box
                        # since it will be different for many boxes
                        connected = connected.astype(float)/len(pt)

                        # fill in all rows of array, one for each box
                        mattemp[i,:] = connected
                        # pdb.set_trace()

                    # add together arrays from different dates
                    mat +=mattemp
                    # pdb.set_trace()
                # divide by the number of different dates combined
                mat /= len(Files)
                np.savez(matfile, mat=mat)

def run_with_times():
    '''
    Combine previously calculated files into month/year files.

    This analysis looks at multiple advection times instead of lumping together.
    '''

    # load in paths for coastline boxes
    d = np.load('calcs/coastpaths.npz', encoding='latin1') # use paths in grid space
    paths = d['pathsg']
    d.close()

    # start indices of drifters that start in each box path
    filename = 'calcs/alongcoastconn/inds-in-coast-paths.npz'
    # indices in pts are referenced to the drifters that start
    # in the shelf transport simulations
    if os.path.exists(filename):
        pts = np.load(filename, encoding='latin1')['pts']
    else:
        # do this once and save
        pts = []
        for i,path in enumerate(paths):
            # which points are inside the regions
            pts.append(find(path.contains_points(np.vstack((xg0.flat, yg0.flat)).T).reshape(xg0.shape)))
        np.savez(filename, pts=pts)

    base = 'calcs/alongcoastconn/conn_in_time/'
    if not os.path.exists(base):
        os.mkdir(base)

    # Loop through along-coast boxes to find which other boxes they are connected to
    # and allow for the 5 crossings kept track of
    # time x boxes x boxes. Using a time for every 4 hours of the 30 total days
    ntimes = int(30*(24/4.))
    Files = glob.glob('calcs/alongcoastconn/20??-??-??T0?.npz')

    for File in Files:
        # print(File)
        matfile = base + File.split('/')[-1]
        if os.path.exists(matfile):
            continue

        d = np.load(File)
        # [# of box paths x # drifters that enter a box x 5 (max # of crosses checked for)]
        inbox = d['inbox'] # time in decimal days when a drifter enters a box path
        # outbox = d['outbox'] # time in decimal days when a drifter exists a box path
        # inds are referenced to the drifters in the shelf transport runs
        inds = d['iinside'] # indices from the original drifters corresponding to in/outbox
        d.close()
        # time x boxes x boxes
        # have to add across time to get results from connected
        mattemp = np.eye(inbox.shape[0])[np.newaxis,:].repeat(ntimes, axis=0)
        # For each box, how many drifters go out from it in to another box?
        for i, path in enumerate(paths):

            pt = pts[i] # indices of drifters that start in this box

            # find indices of where indices in pt can be found in inds/inbox/outbox
            # pt is the same as inds[code]
            code = []
            for p in pt:
                code.extend(find(p==inds))

            # grab all times for when a drifter is entering these boxes.
            # ii: indices of coast boxes drifters are interacting with
            # j: indices of drifters
            ii, j = np.where(~np.isnan(inbox[:,code,0]))
            times = inbox[ii, [code[jj] for jj in j], 0]
            # import pdb; pdb.set_trace()

            # loop through coastal boxes that are used here to save drifter arrivals in time
            # they are repeated in ii, so take set to get unique indices
            for ibox in list(set(ii)):
                iibox = ibox==ii  # inds in ii for box identified by ibox (not unique)

                # Get signal: count drifter entering in time for the coast box ibox
                ndrifters, bins = np.histogram(times[iibox], bins=ntimes, range=(0,30))
                bins = bins[1:]  # skip first entry since these are edges of bins

                # fill in rows of array, one for each box
                mattemp[:,i,ibox] = ndrifters

            # Normalize by the number of drifters that started in box path
            mattemp[:,i,:] /= len(pt)

        np.savez(matfile, mat=mattemp, t=bins)


if __name__ == "__main__":
    calc_bayconn()
    # run()
    # run_with_times()
