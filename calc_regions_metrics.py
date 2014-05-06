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

    # ----- Dispersion ----- #

    # --- horizontal diffusivity --- #
    # doturb=0
    File = '../horizontal_diffusivity/tracks/doturb0_ah0/D2overall.npz'
    d = np.load(File)
    Dhordiff_doturb0 = d['D2']; t = d['t']
    dayshordiff = (t-t[0])/(3600.*24)
    d.close()

    # doturb=2
    File = '../horizontal_diffusivity/tracks/doturb2_ah5/D2overall.npz'
    d = np.load(File)
    Dhordiff_doturb2 = d['D2']
    d.close()

    # --- summer location --- #

    ## summer location - part - 30 -- Jan-Feb -- no diffusion
    #Files = glob('calcs/dispersion/summer/part/30/*-0[1-2]-*.npz')
    #D2 = np.zeros(901); nnans = np.zeros(901)
    #for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    t = d['t']
    #    d.close()
    ## average values for this combination
    #Dsummerpart30_0102 = D2/nnans
    #days = (t-t[0])/(3600.*24)
    #np.savez('calcs/dispersion/summer/part/30/dispersion_0102.npz', days=days, D2=Dsummerpart30_0102)
    d = np.load('calcs/dispersion/summer/part/30/dispersion_0102.npz')
    days = d['days']; Dsummerpart30_0102 = d['D2']

    # # summer location - full - 10 -- Jan-Feb -- no diffusion
    # Files = glob('calcs/dispersion/summer/full/10/*-0[1-2]-*.npz')
    # D2 = np.zeros(901); nnans = np.zeros(901)
    # for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    t = d['t']
    #    d.close()
    # # average values for this combination
    # Dsummerfull10_0102 = D2/nnans
    # days = (t-t[0])/(3600.*24)
    # np.savez('calcs/dispersion/summer/full/10/dispersion_0102.npz', days=days, D2=Dsummerfull10_0102)
    d = np.load('calcs/dispersion/summer/full/10/dispersion_0102.npz')
    days = d['days']; Dsummerfull10_0102 = d['D2']

    # # summer location - part - 30 -- Jan-Feb -- with diffusion
    # Files = glob('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/*-0[1-2]-*.npz')
    # D2 = np.zeros(901); nnans = np.zeros(901)
    # for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    t = d['t']
    #    d.close()
    # # average values for this combination
    # Dsummerpart30_0102_doturb2 = D2/nnans
    # days = (t-t[0])/(3600.*24)
    # np.savez('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/dispersion_0102.npz', days=days, D2=Dsummerpart30_0102_doturb2)
    d = np.load('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/dispersion_0102.npz')
    days = d['days']; Dsummerpart30_0102_doturb2 = d['D2']


    ## summer location - part - 30 -- July-Aug -- no diffusion
    #Files = glob('calcs/dispersion/summer/part/30/*-0[7-8]-*.npz')
    #D2 = np.zeros(901); nnans = np.zeros(901)
    #for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    d.close()
    ## average values for this combination
    #Dsummerpart30_0708 = D2/nnans
    #np.savez('calcs/dispersion/summer/part/30/dispersion_0708.npz', days=days, D2=Dsummerpart30_0708)
    d = np.load('calcs/dispersion/summer/part/30/dispersion_0708.npz')
    Dsummerpart30_0708 = d['D2']

    # # summer location - full - 10 -- July-Aug -- no diffusion
    # Files = glob('calcs/dispersion/summer/full/10/*-0[7-8]-*.npz')
    # D2 = np.zeros(901); nnans = np.zeros(901)
    # for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    d.close()
    # # average values for this combination
    # Dsummerfull10_0708 = D2/nnans
    # np.savez('calcs/dispersion/summer/full/10/dispersion_0708.npz', days=days, D2=Dsummerfull10_0708)
    d = np.load('calcs/dispersion/summer/full/10/dispersion_0708.npz')
    Dsummerfull10_0708 = d['D2']

    # # summer location - part - 30 -- July-Aug -- with diffusion
    # Files = glob('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/*-0[7-8]-*.npz')
    # D2 = np.zeros(901); nnans = np.zeros(901)
    # for File in Files:
    #    # sum all values for this combination
    #    d = np.load(File)
    #    D2 += d['D2']
    #    nnans += d['nnans']
    #    d.close()
    # # average values for this combination
    # Dsummerpart30_0708_doturb2 = D2/nnans
    # np.savez('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/dispersion_0708.npz', days=days, D2=Dsummerpart30_0708_doturb2)
    d = np.load('../shelf_transport_doturb2_ah5/calcs/dispersion/summer/part/30/dispersion_0708.npz')
    Dsummerpart30_0708_doturb2 = d['D2']


    # ----- FSLE ----- #

    alpha = np.sqrt(2)

    # distances increasing with factor alpha
    Rs = np.asarray([np.array([0.7])*alpha**i for i in np.arange(20)]).T # in km

    # --- model output from horizontal_diffusivity: doturb=0 --- #
    # doturb=0, ah=0 -- with all pairs not just close ones
    tSave = np.zeros((1,20))
    nnans = np.zeros((1,20))
    Files = glob('../horizontal_diffusivity/tracks/doturb0_ah0/*fsle.npz')
    for File in Files:
        d = np.load(File)
        tSavetemp = d['tSave']
        ind = ~np.isnan(tSavetemp)
        tSave[ind] += tSavetemp[ind]
        nnans[ind] += d['nnans'][ind]
    #    tSave += d['tSave']
    #    nnans += d['nnans']
        d.close()
    lhordiff_doturb0 = 1/((tSave/nnans))
    #ax.loglog(Rs, l.T, '-', lw=4, color='0.6', ms=10)

    # doturb=2, ah=5
    tSave = np.zeros((1,20))
    nnans = np.zeros((1,20))
    Files = glob('../horizontal_diffusivity/tracks/doturb2_ah5/fsle_allpairs/*fsle.npz')
    for File in Files:
        d = np.load(File)
        tSavetemp = d['tSave']
        ind = ~np.isnan(tSavetemp)
        tSave[ind] += tSavetemp[ind]
        nnans[ind] += d['nnans'][ind]
    #    tSave += d['tSave']
    #    nnans += d['nnans']
        d.close()
    lhordiff_doturb2 = 1/((tSave/nnans))


    # --- summer location --- #

    ## summer location - part - 30 -- 0102 -- no diffusion
    #Files = glob('calcs/fsle/summer/part/30/*-0[1-2]-*.npz')
    #tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    #for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    d.close()
    #lsummerpart30_0102 = 1/((tSave/nnans))
    #np.savez('calcs/fsle/summer/part/30/fsle_0102.npz', Rs=Rs, l=lsummerpart30_0102)
    d = np.load('calcs/fsle/summer/part/30/fsle_0102.npz')
    Rs = d['Rs']; lsummerpart30_0102 = d['l']

    # # summer location - part - 30 -- 0102 -- with diffusion
    # Files = glob('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/*-0[1-2]-*.npz')
    # tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    # for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    pdb.set_trace()
    #    d.close()
    # lsummerpart30_0102_doturb2 = 1/((tSave/nnans))
    # pdb.set_trace()
    # np.savez('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/fsle_0102.npz', Rs=Rs, l=lsummerpart30_0102_doturb2)
    ## d = np.load('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/fsle_0102.npz')
    ## Rs = d['Rs']; lsummerpart30_0102_doturb2 = d['l']

    # # summer location - full - 10 -- 0102 -- no diffusion
    # Files = glob('calcs/fsle/summer/full/10/*-0[1-2]-*.npz')
    # tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    # for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    d.close()
    # lsummerfull10_0102 = 1/((tSave/nnans))
    # np.savez('calcs/fsle/summer/full/10/fsle_0102.npz', Rs=Rs, l=lsummerfull10_0102)
    d = np.load('calcs/fsle/summer/full/10/fsle_0102.npz')
    Rs = d['Rs']; lsummerfull10_0102 = d['l']

    ## summer location - part - 30 -- 0708 -- no diffusion
    #Files = glob('calcs/fsle/summer/part/30/*-0[7-8]-*.npz')
    #tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    #for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    d.close()
    #lsummerpart30_0708 = 1/((tSave/nnans))
    #np.savez('calcs/fsle/summer/part/30/fsle_0708.npz', Rs=Rs, l=lsummerpart30_0708)
    d = np.load('calcs/fsle/summer/part/30/fsle_0708.npz')
    lsummerpart30_0708 = d['l']

    # # summer location - part - 30 -- 0708 -- with diffusion
    # Files = glob('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/*-0[7-8]-*.npz')
    # tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    # for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    d.close()
    # lsummerpart30_0708_doturb2 = 1/((tSave/nnans))
    # np.savez('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/fsle_0708.npz', Rs=Rs, l=lsummerpart30_0708_doturb2)
    d = np.load('../shelf_transport_doturb2_ah5/calcs/fsle/summer/part/30/fsle_0708.npz')
    lsummerpart30_0708_doturb2 = d['l']

    # # summer location - full - 10 -- 0708 -- no diffusion
    # Files = glob('calcs/fsle/summer/full/10/*-0[7-8]-*.npz')
    # tSave = np.zeros((1,20)); nnans = np.zeros((1,20))
    # for File in Files:
    #    d = np.load(File)
    #    tSavetemp = d['tSave']
    #    ind = ~np.isnan(tSavetemp)
    #    tSave[ind] += tSavetemp[ind]
    #    nnans[ind] += d['nnans'][ind]
    #    d.close()
    # lsummerfull10_0708 = 1/((tSave/nnans))
    # np.savez('calcs/fsle/summer/full/10/fsle_0708.npz', Rs=Rs, l=lsummerfull10_0708)
    d = np.load('calcs/fsle/summer/full/10/fsle_0708.npz')
    lsummerfull10_0708 = d['l']


    ## ----- Make the plot ----- ##
    D = np.loadtxt('lacasce_dispersion_Points.txt')

    fig = plt.figure(figsize=(8,12))

    ## Dispersion ##

    ax = fig.add_subplot(2,1,1)
    ax.semilogy(days, Dsummerpart30_0102, '-', color='0.2', lw=3)
    # ax.semilogy(days, Dsummerfull10_0102, ':', color='0.2', lw=3)
    # ax.semilogy(days, Dsummerpart30_0102_doturb2, '--', color='0.2', lw=3)
    ax.semilogy(days, Dsummerpart30_0708, '-', color='0.4', lw=3)
    # ax.semilogy(days, Dsummerfull10_0708, ':', color='0.4', lw=3)
    # ax.semilogy(days, Dsummerpart30_0708_doturb2, '--', color='0.4', lw=3)

    ax.semilogy(dayshordiff, Dhordiff_doturb0, '-.', color='r', lw=3, alpha=0.5)
    ax.semilogy(dayshordiff, Dhordiff_doturb2, '-.', color='r', lw=3, alpha=0.5)

    tind2 = find(days>=2)[0]
    tind9 = find(days>=9)[0]
    tind10 = find(days>=10)[0]
    tind11 = find(days>=11)[0]
    tind24 = find(days>=24)[0]
    ax.semilogy(days[tind2:tind9], .65*np.exp(days[tind2:tind9]*0.55), ':', color='r', lw=2, alpha=0.5)
    ax.semilogy(days[tind11:tind24], 1*days[tind11:tind24]**2.2, ':', color='r', lw=2, alpha=0.5)
    ax.semilogy(D[:,0], D[:,1], '*', color='r', ms=12)

    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Mean separation distance [km$^2\!$]')

    ax.text(12, 4.0, 'Summer', color='0.4')
    ax.text(12, 2.5, 'Winter', color='0.2')
    ax.text(12, 1.5, 'Data/model comparison', color='r')

    ax.set_xlim(0,25)
    ax.set_ylim(1,2e4)
    ax.set_frame_on(False)

    ## FSLE ##

    l = np.loadtxt('LaCasce2008_fsle.txt')

    ax2 = fig.add_subplot(2,1,2)
    ind = ~np.isnan(lsummerpart30_0102)
    ax2.loglog(Rs[ind], lsummerpart30_0102[ind], '-', color='0.2', lw=3)
    # ind = ~np.isnan(lsummerfull10_0102)
    # ax2.loglog(Rs[ind], lsummerfull10_0102[ind], ':', color='0.2', lw=3)
    # ind = ~np.isnan(lsummerpart30_0102_doturb2)
    # ax2.loglog(Rs[ind], lsummerpart30_0102_doturb2[ind], '--', color='0.2', lw=3)

    ind = ~np.isnan(lsummerpart30_0708)
    ax2.loglog(Rs[ind], lsummerpart30_0708[ind], '-', color='0.4', lw=3)
    # ind = ~np.isnan(lsummerfull10_0708)
    # ax2.loglog(Rs[ind], lsummerfull10_0708[ind], ':', color='0.4', lw=3)
    # ind = ~np.isnan(lsummerpart30_0708_doturb2)
    # ax2.loglog(Rs[ind], lsummerpart30_0708_doturb2[ind], '--', color='0.4', lw=3)

    ind = ~np.isnan(lhordiff_doturb0)
    ax2.loglog(Rs[ind], lhordiff_doturb0[ind], '-.', color='r', lw=3, alpha=0.5)
    ind = ~np.isnan(lhordiff_doturb2)
    ax2.loglog(Rs[ind], lhordiff_doturb2[ind], '-.', color='r', lw=3, alpha=0.5)
    #pdb.set_trace()
    ax2.loglog(l[:,0], l[:,1], '*', color='r', ms=12)
    #ind = ~np.isnand(
    ax2.loglog(Rs[0,-11:-4], 1*(Rs[0,-11:-4])**(-2/3.), 'r:', alpha=0.5, lw=2)

    ax2.set_xlabel('Distance [km]')
    ax2.set_ylabel(r'FSLE: $\langle T \rangle ^{-1}$ [days$^{-1}$]')
    ax2.set_frame_on(False)
    ax2.set_xlim(0.6,)

    plt.show()

    # pdb.set_trace()

    fig.savefig('figures/D2_fsle.pdf', bbox_inches='tight', dpi=100)#300)
    fig.savefig('figures/D2_fsle.png', bbox_inches='tight', dpi=300)#300)


def run_dispersion():
    '''
    Run dispersion calculation several areas for shelf transport drifter simulations.
    '''

    # area: 'winter' (Texas) or 'summer' (west LA)
    area = 'summer'

    # amount of area: 'full' or 'part'
    amount = 'full'

    # contour amount: '10', '20', or '30' percent
    contour = '10'

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
    amount = 'full'

    # contour amount: '10', '20', or '30' percent
    contour = '10'

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

                if os.path.exists(fname):
                    continue

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
                ddrifter = 1 # how many drifter indices to include at once
                driftercount = 0

                # logic for looping through more than 1 drifter at once
                while driftercount < ndrifters:
                #    print 'drifter ' + str(driftercount) + ' of ' + str(ndrifters)

                    tstart = time.time()
                    tSavetemp = tracpy.calcs.calc_fsle(lonp[driftercount:driftercount+ddrifter,:], 
                                        latp[driftercount:driftercount+ddrifter,:], tp)
                    ind = ~np.isnan(tSavetemp)
                    tSave += np.nansum(tSavetemp, axis=0)
                    nnans += ind.sum(axis=0)
                    driftercount += ddrifter
                 #   print 'delta time ', time.time()-tstart

                # Save fsle for each file/area combination, NOT averaged
                np.savez(fname, tSave=tSave, nnans=nnans)
                print 'saved file ', fname



if __name__ == "__main__":
    run_fsle()    
