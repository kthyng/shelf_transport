'''
Calculate transport using histograms from make_plots.py.
'''

import numpy as np
import matplotlib.pyplot as plt
import tracpy
import op
import matplotlib as mpl

mpl.rcParams.update({'font.size': 18})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

shelf_depth = 100

whichtime = 'interannual' # 'seasonal' or 'interannual'
whichseason = 'summer' # 'summer' or 'winter'

# Find depths at the center of each histogram bin
loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)


if whichtime == 'seasonal':
    
    # Load in histogram
    d = np.load('figures/cross/seasonal100H.npz')
    H = d['H']
    X, Y = np.meshgrid(op.resize(d['xe'],0), op.resize(d['ye'],0))
    fh = grid['trir'].nn_interpolator(grid['h'].flatten())
    depths = fh(X,Y)

    ishallow = depths < shelf_depth
    ideep = depths > shelf_depth

    # winter
    offshore = np.nansum(H[0,ishallow])/ishallow.sum() # likelihood per histogram bin
    onshore = np.nansum(H[0,ideep])/ideep.sum()


elif whichtime == 'interannual':

    if whichseason == 'winter':

        # Load in histogram
        d = np.load('figures/cross/interannual-winter100H.npz')
        H = d['H']
        X, Y = np.meshgrid(op.resize(d['xe'],0), op.resize(d['ye'],0))
        fh = grid['trir'].nn_interpolator(grid['h'].flatten())
        depths = fh(X,Y)
        flonr = grid['trir'].nn_interpolator(grid['lonr'].flatten())
        lons = flonr(X,Y)
        flatr = grid['trir'].nn_interpolator(grid['latr'].flatten())
        lats = flatr(X,Y)
        d.close()

        ishallow = depths < shelf_depth
        ideep = depths > shelf_depth

        # Calculate the transport for each year for the winter wind transport area.
        # Want places where the depth is less than 100 meters and west of 90 deg.
        iwind = ishallow * (lons<-90) * (lats>27.5)
        Hwind = np.nansum(H[:,iwind.T], axis=1)

        # Plot against calculated mean wind direction
        years = np.arange(2004,2011)
        dirs = np.empty(years.size)
        for i,year in enumerate(years):
            d = np.load('../txla_plots/calcs/wind_stress/calcs' + str(year) + 'winter.npz')
            dirs[i] = d['angle']
            d.close()

        Hwind_rel = (Hwind - Hwind.mean())/(Hwind.max() - Hwind.min())
        fig = plt.figure(figsize=(8.1375,6.9))
        ax = fig.add_subplot(111)
        ax.plot(Hwind_rel, dirs, 'o', color='0.2', ms=10)
        ax.set_xlabel('Transport relative to mean [%]')
        ax.set_ylabel('Mean winter wind direction [deg]')
        fig.savefig('figures/cross/interannual-winter-transport-vs-wind.pdf', bbox_inches='tight')
        np.savez('calcs/shelfconn/' + whichtime + '-' + whichseason + 'transportsum.npz', transport=Hwind_rel)

    elif whichseason == 'summer':

        # Load in histogram
        d = np.load('figures/cross/interannual-summer100H.npz')
        H = d['H']
        X, Y = np.meshgrid(op.resize(d['xe'],0), op.resize(d['ye'],0))
        fh = grid['trir'].nn_interpolator(grid['h'].flatten())
        depths = fh(X,Y)
        flonr = grid['trir'].nn_interpolator(grid['lonr'].flatten())
        lons = flonr(X,Y)
        flatr = grid['trir'].nn_interpolator(grid['latr'].flatten())
        lats = flatr(X,Y)
        d.close()

        ishallow = depths < shelf_depth
        ideep = depths > shelf_depth

        # Calculate the transport for each year for the winter wind transport area.
        # Want places where the depth is less than 100 meters and west of 90 deg.
        iwind = ishallow * (lons<-90) * (lats>27.5)
        Hwind = np.nansum(H[:,iwind.T], axis=1)
        Hwind_rel = (Hwind - Hwind.mean())/(Hwind.max() - Hwind.min())
        np.savez('calcs/shelfconn/' + whichtime + '-' + whichseason + 'transportsum.npz', transport=Hwind_rel)

        # plot vs. wind direction and discharge from Forest et al paper
        # river discharge, m^3/s
        Q = np.array([25800.,17500.,19500.,25400.,44400.,42000.,31300.])
        # wind speed, m/s
        U = np.array([-0.141,-0.472,-1.614,0.426,0.781,2.054,-2.240])
        V = np.array([2.30,0.67,1.81,1.29,1.48,0.91,1.66])
        Dir = np.rad2deg(np.arctan2(V,U))
        S = np.sqrt(U**2+V**2)

        # wind angle and discharge vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Q, Dir, c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.

        # wind magnitude and discharge vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Q, S, c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.

        # wind magnitude and angle vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Dir, S, c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.

        # U wind and discharge vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Q, U, c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.

        # V wind and discharge vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Q, V, c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.

        # U wind * cos(angle) and discharge vs. transport magnitude
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(Q, U*np.cosd(Dir), c=Hwind_rel, cmap='Reds', s=300)#, linewidths=0.
