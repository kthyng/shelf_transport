'''
Calculate the transport from different transport regions in time.
'''

import numpy as np
from glob import glob
import tracpy


ciso = 100  # critical isobath
iciso = 2  # critical isobath index in analysis
region = 3  # Which transport region

Files = glob('calcs/shelfconn/*gc.npz')


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)

d = np.load('calcs/xyp0.npz')
xp0 = d['xp0']; yp0 = d['yp0']
lonp0, latp0 = grid['basemap'](xp0, yp0, inverse=True)
# X, Y = np.meshgrid(op.resize(d['xe'],0), op.resize(d['ye'],0))
fh = grid['trir'].nn_interpolator(grid['h'].flatten())
depths = fh(xp0, yp0)
# flonr = grid['trir'].nn_interpolator(grid['lonr'].flatten())
# lons = flonr(X,Y)
# flatr = grid['trir'].nn_interpolator(grid['latr'].flatten())
# lats = flatr(X,Y)
d.close()

ishallow = depths < ciso
ideep = depths > ciso

# Indices to select out drifters based on where they started in the domain
if region == 1: # shelf bend
    ind = ideep * (lonp0<-95) * (lonp0>-97) * (latp0>26) * (latp0<28)
elif region == 2: # Mississippi region
    ind = ishallow * (lonp0<-88) * (lonp0>-89.3) * (latp0>28.7) * (latp0<30)
elif region == 3: # winter wind transport region
    ind = ishallow * (lonp0<-90) * (latp0>27.5)
elif region == 4: # summer river transport region
    ind = ishallow * (lonp0<-90) * (latp0>27.5)

# ind: drifters that cross the critical isobath that started in the desired region.
# cross: drifters that started in desired region and looking at critical isobath ciso
cross = np.load(Files[0])['cross'][iciso,:][ind]
# innan: indices of drifters that cross ciso
innan = ~np.isnan(cross)
# cross[innan] is drifters that start from the desired region and cross ciso
# Hold time in days that the drifter crossed.
hist, bin_edges = np.histogram(cross[innan], bins=60, range=(0, 30))
