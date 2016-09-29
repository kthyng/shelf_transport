'''
Read in connectivity calculations from find_coastal_path_connectivity.py
and make plots.
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
# import op
from matplotlib.mlab import find
from matplotlib import ticker, colors, cbook
import calendar
import cmocean.cm as cmo

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


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
proj = tracpy.tools.make_proj('nwgom')
grid = tracpy.inout.readgrid(loc, proj)

# load in paths for coastline boxes
d = np.load('calcs/coastpaths.npz') # use paths in grid space
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

X, _ = np.meshgrid(dist, dist)

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
    for i, path in enumerate(pathsxy):
        # only plot every 50th
        if np.mod(i,50)==0:
            if i<110:
                dx = 11000; dy = 0; # shift
            elif i==150:
                dx = 5000; dy = -30000
            elif i>150 and i<300:
                dx = 0; dy = -50000
            elif i==300:
                dy = -40000
            ax.text(path.vertices[0,0] + dx, path.vertices[0,1] + dy, str(i), fontsize=15, color='r')

    fig.savefig('figures/alongcoastconn/domain.png', bbox_inches='tight')


def plot_seasonal():
    '''
    Use calculated files from run() to plot connectivity matrix.
    '''

    cmap = cmo.speed  # 'YlGn'
    log = False

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

    ## Plot setup ##
    fig, axarr = plt.subplots(1,2)#, sharex=True)
    fig.set_size_inches(13, 6.6125)
    fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)
    for i, ax in enumerate(axarr):
       # Titles for subplots
        if i==0:
            # tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title('Winter')
            ax.set_ylabel('Source box number')
            ax.set_xlabel('Destination box number')
            ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
            ax.yaxis.set_ticks_position('left')
            # ax.axis('equal')
            # plt.xticks(xticklocs);
            # plt.yticks(xticklocs);
            # ax.pcolormesh(mat[i,:,:], cmap=cmap)
        elif i==1:
            # tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            ax.set_xlabel('Destination box number')
            ax.xaxis.set_ticks_position('bottom')  # turns off top tick marks
            ax.yaxis.set_ticks_position('left')
            # ax.axis('equal')
            # ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            # plt.xticks(xticklocs);
            # plt.yticks(xticklocs);
            # ax.get_xaxis().set_visible(False)
        if log:
            mappable = ax.pcolormesh(mat[i,:,:]*100., cmap=cmap, vmin=1., vmax=100., norm=colors.LogNorm())
        else:
            mappable = ax.pcolormesh(X, X, mat[i,:,:]*100., cmap=cmap, vmax=100.)
        import pdb; pdb.set_trace()
        ax.set_frame_on(False)
    cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
    if log:
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal', extend='min')
        cb.set_label('Connectivity [%]')
        fig.savefig('figures/alongcoastconn/seasonal-log.png', bbox_inches='tight')
    else:
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')#, pad=0.18)
        cb.set_label('Connectivity [%]')
        fig.savefig('figures/alongcoastconn/seasonal.png', bbox_inches='tight')
    ####





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
    mat = np.zeros((len(paths),len(paths)))
    years = np.arange(2014,2015)
    months = [1,2,7,8]
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

#
# if __name__ == "__main__":
#     run()
