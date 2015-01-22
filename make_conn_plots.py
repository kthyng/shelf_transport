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
import op
from matplotlib.mlab import find
from matplotlib import ticker
from matplotlib import cbook
import calendar

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


def run():

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
    years = np.arange(2004,2015)
    months = [1,2,7,8]
    for year in years:
        for month in months:
            matfile = 'calcs/alongcoastconn/conn-' + str(year) + '-' + str(month).zfill(2) + '.npz'
            if not os.path.exists(matfile):
                Files = glob.glob('calcs/alongcoastconn/' + str(year) \
                            + '-' + str(month).zfill(2) + '-*T00.npz')

                for File in Files:
                    print File
                    d = np.load(File)
                    # [# of box paths x # drifters that enter a box x 5 (max # of crosses checked for)]
                    inbox = d['inbox'] # time in decimal days when a drifter enters a box path
                    outbox = d['outbox'] # time in decimal days when a drifter exists a box path
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

                        pt = pts[i]

                        # find indices of where indices in pt can be found in inds/inbox/outbox
                        # pt is the same as inds[code]
                        code = []
                        for p in pt:
                            code.extend(find(p==inds))

                        # which paths are connected to the box in question
                        connected = ~np.isnan(np.nansum(inbox[:,code,0], axis=1))

                        # put a 1 for the boxes where a drifter will connect
                        mattemp[i,connected] = 1

                    mat +=mattemp
                    # pdb.set_trace()
                np.savez(matfile, mat=mat)


if __name__ == "__main__":
    run()     
