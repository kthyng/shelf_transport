'''
Calculate correlations between alongcoast connectivity (seasonal, interannual)
and forcing mechanisms along the coast.
'''

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import tracpy
from matplotlib.path import Path
import cartopy
import os
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import pandas as pd


def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


## read in alongcoast connectivity line (summed across matrix) ##
# want seasonal and interannual
# work from jupyter notebook in shelf_transport/notebooks/porta_tofrom_surfside.ipynb
base = 'calcs/alongcoastconn/conn_in_time/pa_ss/'
t = np.load(base + 't.npz')['t']  # in decimal days
oname = base + 'aggregated.npz'
Files = glob(base + '2*.npz')
dates = [pd.Timestamp(File.split('/')[-1][:-4]) for File in Files]
d = np.load(oname)
# to save results, [start datetime (dates), 30 days every 4 hours]
ss2pa = d['ss2pa']
pa2ss = d['pa2ss']
d.close()


## read in alongcoast geometric vectors ##
# from shelf_transport/calc_alongcoast_windangle.py
vecname = 'calcs/along_coast_wind/coast_vectors.npz'
d = np.load(vecname)
# upcoast-pointing vector, goes from MX to LA. + is upcoast and toward land.
veccoast = d['veccoast']
# vector that is rotated 90 degrees off from coast direction
veccoast90 = d['veccoast90']
# along-coast distance (km)
dist = d['dist']
# xgcoast=xgcoast, ygcoast=ygcoast, loncoast=loncoast, latcoast=latcoast)


## load in grid
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid('../grid.nc', proj)
# lon_rho = grid.lon_rho[1:-1,1:-1]
# lat_rho = grid.lat_rho[1:-1,1:-1]
# x_rho, y_rho = proj(lon_rho, lat_rho)
# angle = grid.angle_rho[1:-1,1:-1]  # theta to rotate wind vectors


## Select model grid points that are within the coastal boxes ##
inds_rho_boxes = np.load('calcs/coastpaths_pts.npz')['inds_rho_boxes']
inds_u_boxes = np.load('calcs/coastpaths_pts.npz')['inds_u_boxes']
inds_v_boxes = np.load('calcs/coastpaths_pts.npz')['inds_v_boxes']


## read in or calculate along and across coast wind component (seasonal and interannual)
# modified from shelf_transport/calc_alongcoast_windangle.py
# also do river water at the same time while looping

base_wind = 'calcs/along_coast_wind/'
base_river = 'calcs/along_coast_rivers/'

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
m = xr.open_dataset(loc)

years = np.arange(2004, 2015)

for year in years:
    dates = m['ocean_time'].sel(ocean_time=str(year))
    fname_wind = base_wind + str(year) + '.csv'
    fname_river = base_river + str(year) + '.csv'

    if (not os.path.exists(fname_wind)) and (not os.path.exists(fname_river)):

        # units of m^2/s for wind stress
        walong = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
        wacross = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))

        miss = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
        atch = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
        braz = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))


        for i, date in enumerate(dates[::4]):

            datepd = pd.Timestamp(date.values)

            print(datepd)

            # u, v are along and across-shore, surface stress along and across coast
            # these read in an array of model values instead of a vector, so take first row
            # u, v are along and across-shore; take mean since different numbers of grid
            # points in boxes for different grids
            u = [m.sustr.sel(ocean_time=date).isel(eta_u=inds_u_box[0], xi_u=inds_u_box[1])[0,:].mean() for inds_u_box in inds_u_boxes]
            v = [m.svstr.sel(ocean_time=date).isel(eta_v=inds_v_box[0], xi_v=inds_v_box[1])[0,:].mean() for inds_v_box in inds_v_boxes]

            # rotate u, v to be east/north/Cartesian so that I can project properly to coast vectors
            # use mean angle from rho grid since sustr and svstr are on different grids
            angle = [m.angle.isel(eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1]).mean() for inds_rho_box in inds_rho_boxes]
            u, v = rot2d(u, v, angle)  # in cartesian coords after this step

            # make vector of wind
            w = np.vstack((u, v))

            # dot product of the two vectors: gives size of along-shore component
            walong.loc[datepd,:] = np.dot(veccoast, w).diagonal()
            wacross.loc[datepd,:] = np.dot(veccoast90, w).diagonal()

            temp = [m.dye_02.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
            miss.loc[datepd,:] = np.vstack(temp)[:,0]

            temp = [m.dye_03.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
            atch.loc[datepd,:] = np.vstack(temp)[:,0]

            temp = [m.dye_04.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
            braz.loc[datepd,:] = np.vstack(temp)[:,0]

        walong.to_csv(fname_wind)
        wacross.to_csv(fname_wind)
        miss.to_csv(fname_river)
        atch.to_csv(fname_river)
        braz.to_csv(fname_river)


# # calculate percent river amount for all 3 rivers along the coastal boxes
# for year in years:
#     dates = m['ocean_time'].sel(ocean_time=str(year))
#     fname = base + str(year) + '.npz'
#
#     if not os.path.exists(fname):
#
#         miss = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
#         atch = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
#         braz = pd.DataFrame(index=dates[::4], columns=np.arange(0,len(dist)))
#
#         for i, date in enumerate(dates[::4]):
#             datepd = pd.Timestamp(date.values)
#             # if np.mod(i,30*24) == 0:  # print every 30 days
#             print(datepd)
#
#             temp = [m.dye_02.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
#             miss.loc[datepd,:] = np.vstack(temp)[:,0]
#
#             temp = [m.dye_03.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
#             atch.loc[datepd,:] = np.vstack(temp)[:,0]
#
#             temp = [m.dye_04.sel(ocean_time=date).isel(s_rho=-1, eta_rho=inds_rho_box[0], xi_rho=inds_rho_box[1])[0,:].mean() for inds_rho_box in inds_rho_boxes]
#             braz.loc[datepd,:] = np.vstack(temp)[:,0]
#
#         # np.savez(fname, walong=walong, wacross=wacross, dist=dist, dates=dates.data)
#         miss.to_csv(fname)
#         atch.to_csv(fname)
#         braz.to_csv(fname)
