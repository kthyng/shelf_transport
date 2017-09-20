'''

'''

import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry
import shapely.ops
import cartopy
import cartopy.io.shapereader as shpreader
import xarray as xr
import cmocean.cm as cmo
import cartopy
ccrs = cartopy.crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import os
from glob import glob
# from pyproj import Proj

merc = ccrs.Mercator(central_longitude=-85.0)
pc = ccrs.PlateCarree()

grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = xr.open_dataset(grid_filename)

# make Line from main 100 m isobath
cs = plt.contour(grid.lon_rho, grid.lat_rho , grid.h, [100])
p = cs.collections[0].get_paths()[0]  # pull out main isobath (not islands)
v = p.vertices
lon = v[:,0]
lat = v[:,1]
Iso = shapely.geometry.LineString(zip(lon, lat))  # isobath as Line


## not using box intersection anymore ##
# # make Polygons from grid cells bounded by psi grid, and see which intersect
# # the 100 m isobath
# fname = 'calcs/transport/intersecting_boxes.npz'
# if os.path.exists(fname):
#     d = np.load(fname)
#     cells = d['cells']; inds = d['inds']
# else:
#     cells = []
#     jpsi, ipsi = grid.lon_psi.shape
#     inds = []
#     for j in range(0, jpsi-1):
#         for i in range(0, ipsi-1):
#             lattemp = np.hstack((grid.lat_psi[j:j+2,i].data.flatten(),
#                                  grid.lat_psi[j:j+2,i+1][::-1].data.flatten()))
#             lontemp = np.hstack((grid.lon_psi[j:j+2,i].data.flatten(),
#                                  grid.lon_psi[j:j+2,i+1][::-1].data.flatten()))
#             psis = list(zip(lontemp, lattemp))
#             poly = shapely.geometry.Polygon(psis)
#             if poly.intersects(Iso):
#                 cells.append(poly)
#                 inds.append((j,i))
#     np.savez(fname, cells=cells, inds=inds, lon_psi=grid.lon_psi,
#              lat_psi=grid.lat_psi, Iso=Iso)

# # plot to check that boxes are intersecting isobath (previously confirmed)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection=merc)
# lon, lat = list(zip(*Iso.coords))
# ax.plot(lon, lat, transform=pc)
# for cell in cells:
#     lon, lat = list(zip(*cell.boundary.coords))
#     ax.plot(lon, lat, transform=pc)
#


# # Pull out where drifter crosses isobath at saved crossing times
#
# # read in histogram of drifter locations of drifters that cross 100m isobath,
# # calculated in plot_crossshelf_transport.py
# whichseason = 'summer'
# years = np.arange(2004, 2015)
# # for year in years:
# year = years[0]
# fname = 'calcs/shelfconn/crossshelftransport-' + whichseason + '-' + str(year) + '.npz'
# d = np.load(fname)
# H = d['H']; lon_psi = d['lon_psi']; lat_psi = d['lat_psi']
# d.close()
#
# # calculate distance along the isobath
# # inputs = {'proj': 'lcc', 'ellps': 'clrk66', 'datum': 'NAD27',
# #           'lat_1': 22.5, 'lat_2': 31.0, 'lat_0': 30, 'lon_0': -94,
# #           'x_0': 0, 'y_0': 0}
# # proj = Proj(**inputs)
# jbox, ibox = zip(*inds)  # indices of psi grid boxes that intersect 100m isobath
# xbox = grid.x_psi.data[jbox,ibox]
# ybox = grid.y_psi.data[jbox,ibox]
# # distance in km along 100m isobathbox
# d = np.zeros(len(inds))
# for i in range(1,len(inds)):
#     d[i] = d[i-1] + np.sqrt((xbox[i] - xbox[i-1])**2 +
#                             (ybox[i] - ybox[i-1])**2 )/1000.
#
# # use boxes in H that correspond to indices in inds
# plt.figure()
# plt.plot(d, H[jbox,ibox])
# # H[j,i]
