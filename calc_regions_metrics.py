'''
Calculate metrics of drifters starting from specific areas.
'''

import numpy as np
import pdb
from matplotlib.mlab import find
import netCDF4 as netCDF
from scipy import ndimage
import time
from glob import glob
import tracpy
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import tracpy.calcs

# mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


def plot():


    # ----- all the lines ----- #

    # iWTX - 2007 - 02
    Files = glob('calcs/2007-02-*D2iWTX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        t = d['t']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200702W = D2/nnans
    days = (t-t[0])/(3600.*24)

    # iWTX - 2007 - 06
    Files = glob('calcs/2007-06-*D2iWTX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        t = d['t']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200706W = D2/nnans
    days = (t-t[0])/(3600.*24)

    # iWTX - 2008 - 02
    Files = glob('calcs/2008-02-*D2iWTX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        t = d['t']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200802W = D2/nnans
    days = (t-t[0])/(3600.*24)

    # iWTX - 2008 - 06
    Files = glob('calcs/2008-06-*D2iWTX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        t = d['t']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200806W = D2/nnans
    days = (t-t[0])/(3600.*24)

    # iETX - 2007 - 02
    Files = glob('calcs/2007-02-*D2iETX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200702E = D2/nnans

    # iETX - 2007 - 06
    Files = glob('calcs/2007-06-*D2iETX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200706E = D2/nnans

    # iETX - 2008 - 02
    Files = glob('calcs/2008-02-*D2iETX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200802E = D2/nnans

    # iETX - 2008 - 06
    Files = glob('calcs/2008-06-*D2iETX.npz')
    D2 = np.zeros(901); nnans = np.zeros(901)
    for File in Files:
        # sum all values for this combination
        d = np.load(File)
        D2 += d['D2']
        nnans += d['nnans']
        d.close()
    # average values for this combination
    # pdb.set_trace()
    l200806E = D2/nnans

    ## ----- Make the plot ----- ##
    l = np.loadtxt('lacasce_dispersion_Points.txt')

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111)
    # ax.semilogy(days, l200702W, '-', color='darkcyan', lw=3)
    ax.semilogy(days, l200802E, '-', color='darkcyan', lw=3)
    # ax.semilogy(days, l200802W, '-', color='darkcyan', lw=3)
    ax.semilogy(days, l200702E, '-', color='darkcyan', lw=3)
    # ax.semilogy(days, l200806W, '-', color='orange', alpha=0.6, lw=3)
    # ax.semilogy(days, l200706W, '-', color='orange', alpha=0.6, lw=3)
    ax.semilogy(days, l200706E, '-', color='orange', lw=3)
    ax.semilogy(days, l200806E, '-', color='orange', lw=3)

    tind10 = find(days>=10)[0]
    tind25 = find(days>=25)[0]
    ax.semilogy(days[:tind10], 8*np.exp(days[:tind10]*0.55), '--', color='0.5', lw=2)
    ax.semilogy(days[tind10:tind25], 11*days[tind10:tind25]**2.2, '--', color='0.5', lw=2)
    ax.semilogy(l[:,0], l[:,1], '*', color='0.3', ms=8)

    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Mean separation distance [km$^2\!$]')

    ax.text(15, 12, 'Summer', color='orange')
    ax.text(15, 6, 'Winter', color='darkcyan')
    ax.text(15, 3, 'Data', color='0.3')

    ax.set_ylim(1,)
    ax.set_frame_on(False)

    fig.savefig('figures/D2.pdf', bbox_inches='tight', dpi=300)


def run_dispersion():
    '''
    Run dispersion calculation several areas for shelf transport drifter simulations.
    '''

    # area: 'winter' (Texas) or 'summer' (west LA)
    area = 'summer'

    # amount of area: 'full' or 'part'
    amount = 'part'

    # contour amount: '10', '20', or '30' percent
    contour = '20'

    dinds = np.load('calcs/indsinpath-' + area + '-' + amount + '-' + contour + 'percent.npz')
    inds = dinds['inds']
    dinds.close()

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, usebasemap=False)

    floc = 'calcs/dispersion/' + area + '/' + amount + '/' + contour

    if not os.path.exists(floc):
        os.makedirs(floc)


    # run through each year and month
    for year in np.arange(2004,2011):
        for month in np.arange(1,13):

            Files = glob('tracks/' + str(year) + '-' + str(month).zfill(2) + '-*gc.nc')

            for File in Files:

                # pdb.set_trace()

                fname = floc + '/' + File.split('/')[-1][:-5] + '.npz'
                # pdb.set_trace()
                d = netCDF.Dataset(File)
                xg = d.variables['xg'][:]
                yg = d.variables['yg'][:]
                tp = d.variables['tp'][:]
                d.close()

                # only want drifters starting in this region
                xg = xg[inds,:]; yg = yg[inds,:]

                nanind = np.isnan(xg) + (xg==-1) + (np.ceil(xg)<=5) + (np.ceil(xg)>=grid['xr'].shape[0]-5) + (np.ceil(yg)<=5)
                lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll') 
                lonp[nanind] = np.nan; latp[nanind] = np.nan
                del(xg,yg) # don't need grid info anymore

                nt = tp.size
                ndrifters = lonp.shape[0]


                ## ---- Find pairs of drifters in each region ---- ##

                # let the index in axis 0 be the drifter id
                ID = np.arange(lonp.shape[0])
                # pdb.set_trace()
                # save pairs to save time since they are always the same
                if not os.path.exists(floc + '/pairs.npz'):

                    dist = np.zeros((lonp.shape[0],lonp.shape[0]))
                    for idrifter in xrange(lonp.shape[0]):
                        print 'getting distances:'
                        print 'drifter ' + str(idrifter) + ' of ' + str(ndrifters)
                        # dist contains all of the distances from other drifters for each drifter
                        dist[idrifter,:] = tracpy.calcs.get_dist(lonp[idrifter,0], lonp[:,0], latp[idrifter,0], latp[:,0])

                    pairs = []
                    for idrifter in xrange(lonp.shape[0]):
                        print 'finding pairs'
                        print 'drifter ' + str(idrifter) + ' of ' + str(ndrifters)
                        ind = find(dist[idrifter,:]<=1)
                        for i in ind:
                            if ID[idrifter] != ID[i]:
                                pairs.append([min(ID[idrifter], ID[i]), 
                                                max(ID[idrifter], ID[i])])

                    pairs_set = set(map(tuple,pairs))
                    pairs = map(list,pairs_set)# now pairs has only unique pairs of drifters
                    pairs.sort() #unnecessary but handy for checking work
                    np.savez(floc + '/pairs.npz', pairs=pairs)
                else:
                    pairs = np.load(floc + '/pairs.npz')['pairs']

                ## ---- Finished finding pairs ---- ##

                # Loop over pairs of drifters from this area/time period and sum the FSLE, 
                # then average at the end

                D2 = np.zeros(nt)
                nnans = np.zeros(nt) # to collect number of non-nans over all drifters for a time
                for j, ipair in enumerate(xrange(len(pairs))):

                    if j==0:
                        print 'calculating dispersion for tracks ' + File.split('/')[-1][:-5]

                    dist = tracpy.calcs.get_dist(lonp[pairs[ipair][0],:], lonp[pairs[ipair][1],:], 
                                latp[pairs[ipair][0],:], latp[pairs[ipair][1],:])
                    nnanstemp = ~np.isnan(dist)
                    D2[nnanstemp] += dist[nnanstemp]**2
                    nnans += nnanstemp

                # Save fsle for each file/area combination, NOT averaged
                np.savez(fname, D2=D2, nnans=nnans, t=tp)
                print 'saved file', fname

def run_fsle():
    '''
    Run FSLE calculation several areas for shelf transport drifter simulations.
    '''

    # area: 'winter' (Texas) or 'summer' (west LA)
    area = 'summer'

    # amount of area: 'full' or 'part'
    amount = 'part'

    # contour amount: '10', '20', or '30' percent
    contour = '20'

    dinds = np.load('calcs/indsinpath-' + area + '-' + amount + '-' + contour + 'percent.npz')
    inds = dinds['inds']
    dinds.close()

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc, usebasemap=False)

    floc = 'calcs/fsle/' + area + '/' + amount + '/' + contour

    if not os.path.exists(floc):
        os.makedirs(floc)


    # run through each year and month
    for year in np.arange(2004,2011):
        for month in np.arange(1,13):

            Files = glob('tracks/' + str(year) + '-' + str(month).zfill(2) + '-*gc.nc')

            for File in Files:

                # pdb.set_trace()

                fname = floc + '/' + File.split('/')[-1][:-5] + '.npz'
                # pdb.set_trace()
                d = netCDF.Dataset(File)
                xg = d.variables['xg'][:]
                yg = d.variables['yg'][:]
                tp = d.variables['tp'][:]
                d.close()

                # only want drifters starting in this region
                xg = xg[inds,:]; yg = yg[inds,:]

                nanind = np.isnan(xg) + (xg==-1) + (np.ceil(xg)<=5) + (np.ceil(xg)>=grid['xr'].shape[0]-5) + (np.ceil(yg)<=5)
                lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll') 
                lonp[nanind] = np.nan; latp[nanind] = np.nan
                del(xg,yg) # don't need grid info anymore

                nt = tp.size
                ndrifters = lonp.shape[0]

                # Loop over pairs of drifters from this area/time period and sum the FSLE, 
                # then average at the end

                ndrifters = lonp.shape[0]
                tSave = np.zeros((1,20))
                nnans = np.zeros((1,20)) # to collect number of non-nans over all drifters for a time
                ddrifter = 500 # how many drifter indices to include at once
                driftercount = 0

                # logic for looping through more than 1 drifter at once
                while driftercount < ndrifters:
                    print 'drifter ' + str(driftercount) + ' of ' + str(ndrifters)
                    tSavetemp = tracpy.calcs.calc_fsle(lonp[driftercount:driftercount+ddrifter,:], 
                                        latp[driftercount:driftercount+ddrifter,:], tp)
                    ind = ~np.isnan(tSavetemp)
                    tSave += np.nansum(tSavetemp, axis=0)
                    nnans += ind.sum(axis=0)
                    driftercount += ddrifter

                # Save fsle for each file/area combination, NOT averaged
                np.savez(fname, tSave=tSave, nnans=nnans)
                print 'saved file', fname



if __name__ == "__main__":
    run_fsle()    
