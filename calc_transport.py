'''
Calculate transport using histograms from make_plots.py.
'''

import numpy as np
import matplotlib.pyplot as plt
import tracpy
import op

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

shelf_depth = 100

whichtime = 'interannual' # 'seasonal' or 'interannual'

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
    plt.plot(Hwind_rel, dirs - 180, 'ko', ms=10)
    plt.xlabel('Transport relative to mean [%]')
    plt.ylabel('Mean winter wind direction [deg, relative to 180]')
    plt.savefig('figures/cross/interannual-winter-transport-vs-wind.png', bbox_inches='tight')