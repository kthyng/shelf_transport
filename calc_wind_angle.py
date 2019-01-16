'''
Calculate angle between mean wind for each season relative to coastline
Run with python2
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

fname = 'calcs/alongcoastconn/alongcoastwindangle.npz'
matrixfname = 'calcs/alongcoastconn/conn-seasonal.npz'

# calculate along-coast angle if needed
if not os.path.exists(fname):

    # load in grid
    proj = tracpy.tools.make_proj('nwgom')
    grid = tracpy.inout.readgrid('../grid.nc', proj)
    lonr = grid.lon_rho[1:-1,1:-1]
    latr = grid.lat_rho[1:-1,1:-1]
    xr, yr = proj(lonr, latr)

    # load in coastline as clicked in by coast boxes
    outerpathxy = np.load('calcs/coastpaths.npz', encoding='latin1')['outerpathxy'].item()
    # # x, y pts for path closest to land
    # xcoast = outerpathxy.vertices[:343,0]; ycoast = outerpathxy.vertices[:343,1]
    # x, y pts for path farther from land
    xcoast = outerpathxy.vertices[343:-2,0]; ycoast = outerpathxy.vertices[343:-2,1]
    # also flipping so go from MX to LA in order
    xcoast = xcoast[::-1]; ycoast = ycoast[::-1]

    # read in average wind from surface salinity file, on psi grid I think
    d = np.load('../txla_plots/calcs/means/salt-Winter.npz')
    sustrw = d['sustr']
    svstrw = d['svstr']
    d.close()
    d = np.load('../txla_plots/calcs/means/salt-Summer.npz')
    sustrs = d['sustr']
    svstrs = d['svstr']
    d.close()


    # interpolate wind information onto outerpathxy vertices
    # Triangulation for grid space to curvilinear space, but 1:-1 in x,y for rho grid
    pts = np.column_stack((xr.flatten(), yr.flatten()))
    tess = Delaunay(pts)
    tri = mtri.Triangulation(xr.flatten(), yr.flatten(), tess.simplices.copy())

    # winter
    fsustrw = mtri.LinearTriInterpolator(tri, sustrw.flatten())
    fsvstrw = mtri.LinearTriInterpolator(tri, svstrw.flatten())
    sustrwcoast = fsustrw(xcoast, ycoast)
    svstrwcoast = fsvstrw(xcoast, ycoast)
    # summer
    fsustrs = mtri.LinearTriInterpolator(tri, sustrs.flatten())
    fsvstrs = mtri.LinearTriInterpolator(tri, svstrs.flatten())
    sustrscoast = fsustrs(xcoast, ycoast)
    svstrscoast = fsvstrs(xcoast, ycoast)


    # Calculate angle between interpolated wind angle and coastline to get:
    # -1 for downcoast parallel, 0 for perpendicual, and +1 for upcoast parallel

    # make upcoast-pointing vector from coast coordinates on either side of given pt
    # goes from MX to LA. + is upcoast and toward land.
    veccoast = np.zeros((xcoast.size, 2))*np.nan
    veccoast[1:-1,0] = (xcoast[2:] - xcoast[0:-2])
    veccoast[1:-1,1] = (ycoast[2:] - ycoast[0:-2])
    # unit vector: normalize
    veccoast /= np.sqrt(veccoast[:,0]**2 + veccoast[:,1]**2)[:,np.newaxis]

    # winter
    # normalize wind vector to unit vector
    sstrw = np.vstack((sustrwcoast, svstrwcoast)).T
    sstrw /= np.sqrt(sstrw[:,0]**2 + sstrw[:,1]**2)[:,np.newaxis]
    # dot product of the two vectors
    anglew = np.dot(veccoast, sstrw.T).diagonal()
    # summer
    # normalize wind vector to unit vector
    sstrs = np.vstack((sustrscoast, svstrscoast)).T
    sstrs /= np.sqrt(sstrs[:,0]**2 + sstrs[:,1]**2)[:,np.newaxis]
    # dot product of the two vectors
    angles = np.dot(veccoast, sstrs.T).diagonal()


    # test in plot
    # coast path points
    plt.plot(xcoast, ycoast)
    # winds along coast path
    plt.quiver(xcoast, ycoast, sustrwcoast, svstrwcoast);
    # angle
    plt.scatter(xcoast, ycoast, c=anglew, s=80, cmap=cmo.balance, vmin=-1, vmax=1)
    plt.colorbar()

    # save
    np.savez(fname,
             xcoast=xcoast, ycoast=ycoast,
             sustrwcoast=sustrwcoast.data, svstrwcoast=svstrwcoast.data,
             sustrscoast=sustrscoast.data, svstrscoast=svstrscoast.data,
             anglew=anglew.data, angles=angles.data)

else:  # read in if already done this

    d = np.load(fname)
    anglew = d['anglew']
    angles = d['angles']
    d.close()



# read in along-coast seasonal connectivity info from make_conn_plots.py
mat = np.load(matrixfname)['mat']

# load in paths for coastline boxes
d = np.load('calcs/coastpaths.npz') # use paths in grid space
pathsxy = d['pathsxy']
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
dmax = dist.max()

# sum to calculate how much of the stuff from a box goes up (positive) or downcoast (negative)
# units are percent of possible total drifters starting in a box
# along: 2 seasons x 342 boxes x 2 directions (down and upcoast)
alongstart = np.zeros((2, 342, 2))*np.nan  # add extra x2 to keep downcoast (-) and upcoast (+) separate
for row in range(mat.shape[1]):
    for season in range(mat.shape[0]):
        alongstart[season,row,0] = -mat[season,row,:row+1].sum()  # downcoast
        alongstart[season,row,1] = mat[season,row,row:].sum()  # upcoast
# alongend is summing where drifters end up
alongend = np.zeros((2, 2, 342))*np.nan  # add extra x2 to keep downcoast (-) and upcoast (+) separate
for col in range(mat.shape[2]):
    for season in range(mat.shape[0]):
        alongend[season,0,col] = -mat[season,col:,col].sum()  # downcoast
        alongend[season,1,col] = mat[season,:col+1,col].sum()  # upcoast


colu = '#218983'  # upcoast color
cold = '#cb6863'  # downcoast color
# fig = plt.figure(figsize=(14,6))
fig, axes = plt.subplots(2, 1, figsize=(14,6), sharex=True)
fig.subplots_adjust(left=0.05, right=0.95, hspace=0.12, top=0.97)
# winter
axes[0].grid('on', linestyle='-', alpha=0.5, linewidth=0.1)
axes[0].plot(dist, alongend[0,0,:], '--', color=cold, lw=2)  # downcoast
axes[0].plot(dist, alongend[0,1,:], '--', color=colu, lw=2)  # upcoast
amax = abs(alongend).max()
# add zero line
axes[0].plot(dist, np.zeros(dist.size), color='0.7', lw=3, alpha=0.5)
# ax.plot(dist, alongstart[0,:,0], color='#cb6863', lw=2)  # downcoast
# ax.plot(dist, alongstart[0,:,1], color='#218983', lw=2)  # upcoast
# amax = max((abs(alongstart).max(),abs(alongend).max()))
axes[0].axis('tight')
axes[0].set_ylim(-amax, amax)
axes[0].text(0.13, 0.9, 'Winter', fontsize=16, transform=axes[0].transAxes)
axes[0].set_ylabel('Along-coast connectivity [%]')
# labels
axes[0].text(0.01, 0.85, 'upcoast', color=colu, rotation=90, transform=axes[0].transAxes, alpha=0.6)
axes[0].text(0.01, 0.4, 'downcoast', color=cold, rotation=90, transform=axes[0].transAxes, alpha=0.6)
# add pearson r
ind = ~np.isnan(anglew)
# sum across both up and downcoast for pearson r
r, p = scipy.stats.pearsonr(abs(alongend[0,:,ind]).sum(axis=1),anglew[ind])
axes[0].text(0.605, 0.83, 'Pearson r=%3.2f, p=%2.4f' % (r,p), transform=axes[0].transAxes, color='0.3')
# add angle
ax2 = axes[0].twinx()
ax2.plot(dist, anglew, ':', color='0.2', lw=1)
ax2.axis('tight')
ax2.set_ylim(-1.1, 1.1)
ax2.set_ylabel('Relative angle')

# summer
axes[1].grid('on', linestyle='-', alpha=0.5, linewidth=0.1)
axes[1].plot(dist, alongend[1,0,:], '--', color=cold, lw=2)  # downcoast
axes[1].plot(dist, alongend[1,1,:], '--', color=colu, lw=2)  # upcoast
amax = abs(alongend).max()
# add zero lines
axes[1].plot(dist, np.zeros(dist.size), color='0.7', lw=3, alpha=0.5)
axes[1].axis('tight')
axes[1].set_ylim(-amax, amax)
axes[1].text(0.13, 0.9, 'Summer', fontsize=16, transform=axes[1].transAxes)
axes[1].set_ylabel('Along-coast connectivity [%]')
# labels
axes[1].text(0.01, 0.85, 'upcoast', color=colu, rotation=90, transform=axes[1].transAxes, alpha=0.6)
axes[1].text(0.01, 0.4, 'downcoast', color=cold, rotation=90, transform=axes[1].transAxes, alpha=0.6)
# add pearson r
ind = ~np.isnan(angles)
r, p = scipy.stats.pearsonr(abs(alongend[1,:,ind]).sum(axis=1),angles[ind])
axes[1].text(0.605, 0.83, 'Pearson r=%3.2f, p=%2.4f' % (r,p), transform=axes[1].transAxes, color='0.3')
# add angle
ax2 = axes[1].twinx()
ax2.plot(dist, angles, ':', color='0.2', lw=1)
ax2.axis('tight')
ax2.set_ylim(-1.1, 1.1)
ax2.set_ylabel('Relative angle')
fig.savefig('figures/alongcoastconn/meanangle.pdf', bbox_inches='tight')


# # look at mean salinity along the coast
# saltw = np.load('../txla_plots/calcs/means/salt-Winter.npz')['salt']
# fsaltw = mtri.LinearTriInterpolator(tri, saltw.flatten())
# saltwcoast = fsaltw(xcoast, ycoast)
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# # winter acconn, relative angle, distance, salinity
# ax.scatter(anglew, alongend[0,:,:].sum(axis=0), s=dist/10., c=saltwcoast, cmap=cmo.haline, alpha=0.5)
# # ax.scatter(dist, alongend[0,:,:].sum(axis=0), c=anglew, s=saltwcoast, cmap=cmo.haline, alpha=0.5)


# Look at distance between coastpath and 100 m isobath

# # get locations of 100 m isobath
# plt.figure(); con = plt.contour(grid.x_rho, grid.y_rho, grid.h,[100])
# p100 = con.allsegs[0][0]
# # plot
# plt.plot(xcoast, ycoast, 'k')
# plt.plot(p100[:,0], p100[:,1], 'r')


# try polar plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.scatter(np.arccos(anglew), np.ones(anglew.size),
           c=alongend[0,:,ind].sum(axis=1)[np.newaxis,:],cmap=cmo.curl_r)
