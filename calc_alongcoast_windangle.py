'''
Calculate angle between wind in time for each season relative to coastline
'''

import numpy as np
import matplotlib.pyplot as plt
import tracpy
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
import cmocean.cm as cmo
import os
import matplotlib as mpl
import scipy.stats
import xarray as xr


mpl.rcParams.update({'font.size': 14})

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr


base = 'calcs/along_coast_wind/'
if not os.path.exists(base):
    os.makedirs(base)

# load in grid
proj = tracpy.tools.make_proj('nwgom-pyproj')
grid = tracpy.inout.readgrid('../../grid.nc', proj)
lon_rho = grid.lon_rho[1:-1,1:-1]
lat_rho = grid.lat_rho[1:-1,1:-1]
x_rho, y_rho = proj(lon_rho, lat_rho)
angle = grid.angle_rho[1:-1,1:-1]  # theta to rotate wind vectors

# load in coastline as clicked in by coast boxes
outerpathxy = np.load('calcs/coastpaths.npz', encoding='latin1')['outerpathg'].item()
# # x, y pts for path closest to land
# xcoast = outerpathxy.vertices[:343,0]; ycoast = outerpathxy.vertices[:343,1]
# x, y pts for path farther from land
xgcoast = outerpathxy.vertices[343:-2,0]; ygcoast = outerpathxy.vertices[343:-2,1]
# also flipping so go from MX to LA in order
xgcoast = xgcoast[::-1]; ygcoast = ygcoast[::-1]

# convert from grid coords to xy
xcoast, ycoast, _ = tracpy.tools.interpolate2d(xgcoast, ygcoast, grid, 'm_ij2xy')

# convert from grid coords to ll to save
loncoast, latcoast, _ = tracpy.tools.interpolate2d(xgcoast, ygcoast, grid, 'm_ij2ll')



# make upcoast-pointing vector from coast coordinates on either side of given pt
# goes from MX to LA. + is upcoast and toward land.
veccoast = np.zeros((xcoast.size, 2))*np.nan
veccoast[1:-1,0] = (xcoast[2:] - xcoast[0:-2])
veccoast[1:-1,1] = (ycoast[2:] - ycoast[0:-2])
# unit vector: normalize
veccoast /= np.sqrt(veccoast[:,0]**2 + veccoast[:,1]**2)[:,np.newaxis]

# calculate vector that is rotated 90 degrees off from coast direction
# cos(90) = 0, sin(90) = 1
veccoast90 = np.vstack((veccoast[:,0]*0 - veccoast[:,1]*1,
                       veccoast[:,0]*1 + veccoast[:,1]*0)).T


# along-coast distance
dtemp = np.sqrt((xcoast[1:] - xcoast[:-1])**2 + (ycoast[1:] - ycoast[:-1])**2)
dist = dtemp.cumsum()
dist /= 1000 # convert to km
dist = np.hstack((0, dist))


# save
vecname = base + 'coast_vectors.npz'
np.savez(vecname, veccoast=veccoast, veccoast90=veccoast90, dist=dist,
         xgcoast=xgcoast, ygcoast=ygcoast, loncoast=loncoast, latcoast=latcoast)


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
m = xr.open_dataset(loc)

pts = np.column_stack((x_rho.flatten(), y_rho.flatten()))
tess = Delaunay(pts)
tri = mtri.Triangulation(x_rho.flatten(), y_rho.flatten(), tess.simplices.copy())




years = np.arange(2004, 2015)

for year in years:
    dates = m['ocean_time'].sel(ocean_time=str(year))
    fname = base + str(year) + '.npz'

    if not os.path.exists(fname):

        walong = np.empty((len(dates[::4]), 342)); wacross = np.empty((len(dates[::4]), 342))

        for i, date in enumerate(dates[::4]):
            # if np.mod(i,30*24) == 0:  # print every 30 days
            print(date)

            # U and V wind speed on rho grid but without ghost cells
            # u, v are along and across-shore
            u = m.Uwind.sel(ocean_time=date).isel(eta_rho=slice(1,-1), xi_rho=slice(1,-1))
            v = m.Vwind.sel(ocean_time=date).isel(eta_rho=slice(1,-1), xi_rho=slice(1,-1))

            # rotate u, v to be east/north/Cartesian so that I can project properly to coast vectors
            u, v = rot2d(u, v, angle)

            # interpolate wind information onto outerpathxy vertices
            fu = mtri.LinearTriInterpolator(tri, u.data.flatten())
            fv = mtri.LinearTriInterpolator(tri, v.data.flatten())
            ucoast = fu(xcoast, ycoast)  # x-direction wind, at locations along coast
            vcoast = fv(xcoast, ycoast)  # y-direction wind, at locations along coast

            # make vector of wind
            w = np.vstack((ucoast, vcoast))

            # dot product of the two vectors: gives size of along-shore component
            walong[i,:] = np.dot(veccoast, w).diagonal()
            wacross[i,:] = np.dot(veccoast90, w).diagonal()
        np.savez(fname, walong=walong, wacross=wacross, dist=dist, dates=dates.data)

# # sum to calculate how much of the stuff from a box goes up (positive) or downcoast (negative)
# # units are percent of possible total drifters starting in a box
# # along: 2 seasons x 342 boxes x 2 directions (down and upcoast)
# alongstart = np.zeros((2, 342, 2))*np.nan  # add extra x2 to keep downcoast (-) and upcoast (+) separate
# for row in range(mat.shape[1]):
#     for season in range(mat.shape[0]):
#         alongstart[season,row,0] = -mat[season,row,:row+1].sum()  # downcoast
#         alongstart[season,row,1] = mat[season,row,row:].sum()  # upcoast
# # alongend is summing where drifters end up
# alongend = np.zeros((2, 2, 342))*np.nan  # add extra x2 to keep downcoast (-) and upcoast (+) separate
# for col in range(mat.shape[2]):
#     for season in range(mat.shape[0]):
#         alongend[season,0,col] = -mat[season,col:,col].sum()  # downcoast
#         alongend[season,1,col] = mat[season,:col+1,col].sum()  # upcoast
#
#
# colu = '#218983'  # upcoast color
# cold = '#cb6863'  # downcoast color
# # fig = plt.figure(figsize=(14,6))
# fig, axes = plt.subplots(2, 1, figsize=(14,6), sharex=True)
# fig.subplots_adjust(left=0.05, right=0.95, hspace=0.12, top=0.97)
# # winter
# axes[0].grid('on', linestyle='-', alpha=0.5, linewidth=0.1)
# axes[0].plot(dist, alongend[0,0,:], '--', color=cold, lw=2)  # downcoast
# axes[0].plot(dist, alongend[0,1,:], '--', color=colu, lw=2)  # upcoast
# amax = abs(alongend).max()
# # add zero line
# axes[0].plot(dist, np.zeros(dist.size), color='0.7', lw=3, alpha=0.5)
# # ax.plot(dist, alongstart[0,:,0], color='#cb6863', lw=2)  # downcoast
# # ax.plot(dist, alongstart[0,:,1], color='#218983', lw=2)  # upcoast
# # amax = max((abs(alongstart).max(),abs(alongend).max()))
# axes[0].axis('tight')
# axes[0].set_ylim(-amax, amax)
# axes[0].text(0.13, 0.9, 'Winter', fontsize=16, transform=axes[0].transAxes)
# axes[0].set_ylabel('Along-coast connectivity [%]')
# # labels
# axes[0].text(0.01, 0.85, 'upcoast', color=colu, rotation=90, transform=axes[0].transAxes, alpha=0.6)
# axes[0].text(0.01, 0.4, 'downcoast', color=cold, rotation=90, transform=axes[0].transAxes, alpha=0.6)
# # add pearson r
# ind = ~np.isnan(anglew)
# # sum across both up and downcoast for pearson r
# r, p = scipy.stats.pearsonr(abs(alongend[0,:,ind]).sum(axis=1),anglew[ind])
# axes[0].text(0.605, 0.83, 'Pearson r=%3.2f, p=%2.4f' % (r,p), transform=axes[0].transAxes, color='0.3')
# # add angle
# ax2 = axes[0].twinx()
# ax2.plot(dist, anglew, ':', color='0.2', lw=1)
# ax2.axis('tight')
# ax2.set_ylim(-1.1, 1.1)
# ax2.set_ylabel('Relative angle')
#
# # summer
# axes[1].grid('on', linestyle='-', alpha=0.5, linewidth=0.1)
# axes[1].plot(dist, alongend[1,0,:], '--', color=cold, lw=2)  # downcoast
# axes[1].plot(dist, alongend[1,1,:], '--', color=colu, lw=2)  # upcoast
# amax = abs(alongend).max()
# # add zero lines
# axes[1].plot(dist, np.zeros(dist.size), color='0.7', lw=3, alpha=0.5)
# axes[1].axis('tight')
# axes[1].set_ylim(-amax, amax)
# axes[1].text(0.13, 0.9, 'Summer', fontsize=16, transform=axes[1].transAxes)
# axes[1].set_ylabel('Along-coast connectivity [%]')
# # labels
# axes[1].text(0.01, 0.85, 'upcoast', color=colu, rotation=90, transform=axes[1].transAxes, alpha=0.6)
# axes[1].text(0.01, 0.4, 'downcoast', color=cold, rotation=90, transform=axes[1].transAxes, alpha=0.6)
# # add pearson r
# ind = ~np.isnan(angles)
# r, p = scipy.stats.pearsonr(abs(alongend[1,:,ind]).sum(axis=1),angles[ind])
# axes[1].text(0.605, 0.83, 'Pearson r=%3.2f, p=%2.4f' % (r,p), transform=axes[1].transAxes, color='0.3')
# # add angle
# ax2 = axes[1].twinx()
# ax2.plot(dist, angles, ':', color='0.2', lw=1)
# ax2.axis('tight')
# ax2.set_ylim(-1.1, 1.1)
# ax2.set_ylabel('Relative angle')
# fig.savefig('figures/alongcoastconn/meanangle.pdf', bbox_inches='tight')
#
#
# # # look at mean salinity along the coast
# # saltw = np.load('../txla_plots/calcs/means/salt-Winter.npz')['salt']
# # fsaltw = mtri.LinearTriInterpolator(tri, saltw.flatten())
# # saltwcoast = fsaltw(xcoast, ycoast)
# #
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # # winter acconn, relative angle, distance, salinity
# # ax.scatter(anglew, alongend[0,:,:].sum(axis=0), s=dist/10., c=saltwcoast, cmap=cmo.haline, alpha=0.5)
# # # ax.scatter(dist, alongend[0,:,:].sum(axis=0), c=anglew, s=saltwcoast, cmap=cmo.haline, alpha=0.5)
#
#
# # Look at distance between coastpath and 100 m isobath
#
# # # get locations of 100 m isobath
# # plt.figure(); con = plt.contour(grid.x_rho, grid.y_rho, grid.h,[100])
# # p100 = con.allsegs[0][0]
# # # plot
# # plt.plot(xcoast, ycoast, 'k')
# # plt.plot(p100[:,0], p100[:,1], 'r')
#
#
# # try polar plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='polar')
# ax.scatter(np.arccos(anglew), np.ones(anglew.size),
#            c=alongend[0,:,ind].sum(axis=1)[np.newaxis,:],cmap=cmo.curl_r)
