'''
Use outerpath in lat/lon of coastal boxes to select grid nodes within boxes.
'''

import numpy as np
import netCDF4 as netCDF
import matplotlib.pyplot as plt

## points for full strip of area near coast ##
outerpathll = np.load('calcs/coastpath_outerll.npz', encoding='latin1')['outerpathll'].item()

grid = netCDF.Dataset('../grid.nc')

pts_u = np.vstack((grid.lon_u.flatten(), grid.lat_u.flatten())).T
pts_v = np.vstack((grid.lon_v.flatten(), grid.lat_v.flatten())).T
pts_rho = np.vstack((grid.lon_rho.flatten(), grid.lat_rho.flatten())).T
pts_psi = np.vstack((grid.lon_psi.flatten(), grid.lat_psi.flatten())).T

inds = outerpathll.contains_points(pts_u)
pts_u_all = pts_u[inds]
inds = outerpathll.contains_points(pts_v)
pts_v_all = pts_v[inds]
inds = outerpathll.contains_points(pts_rho)
pts_rho_all = pts_rho[inds]
inds = outerpathll.contains_points(pts_psi)
pts_psi_all = pts_psi[inds]

## points by box ##
pathsll = np.load('calcs/coastpaths.npz', encoding='latin1')['paths']

pts_u_boxes = [pts_u[pathll.contains_points(pts_u)] for pathll in pathsll]
pts_v_boxes = [pts_v[pathll.contains_points(pts_v)] for pathll in pathsll]
pts_rho_boxes = [pts_rho[pathll.contains_points(pts_rho)] for pathll in pathsll]
pts_psi_boxes = [pts_psi[pathll.contains_points(pts_psi)] for pathll in pathsll]

## indices by box ##
# find points in each box, find flattened index, unravel to coordinate (2d) index
inds_u_boxes = [np.unravel_index(np.where(pathll.contains_points(pts_u))[0], grid.lon_u.shape) for pathll in pathsll]
inds_v_boxes = [np.unravel_index(np.where(pathll.contains_points(pts_v))[0], grid.lon_v.shape) for pathll in pathsll]
inds_rho_boxes = [np.unravel_index(np.where(pathll.contains_points(pts_rho))[0], grid.lon_rho.shape) for pathll in pathsll]
inds_psi_boxes = [np.unravel_index(np.where(pathll.contains_points(pts_psi))[0], grid.lon_psi.shape) for pathll in pathsll]

# check - looks good!
plt.plot(*outerpathll.vertices.T)
plt.plot(*pts_u_all.T, 'o');
[plt.plot(*pathll.vertices.T) for pathll in pathsll];
plt.plot(*pts_u_boxes[100].T, 'x');
plt.plot(grid.lon_u[inds_u_boxes[100]], grid.lat_u[inds_u_boxes[100]], '.k');

## save ##
# lon/lat points
np.savez('calcs/coastpaths_pts.npz', pts_u_all=pts_u_all, pts_v_all=pts_v_all,
         pts_rho_all=pts_rho_all, pts_psi_all=pts_psi_all,
         pts_u_boxes=pts_u_boxes, pts_v_boxes=pts_v_boxes,
         pts_rho_boxes=pts_rho_boxes, pts_psi_boxes=pts_psi_boxes,
         inds_u_boxes=inds_u_boxes, inds_v_boxes=inds_v_boxes,
         inds_rho_boxes=inds_rho_boxes, inds_psi_boxes=inds_psi_boxes)
