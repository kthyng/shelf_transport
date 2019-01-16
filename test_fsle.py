'''
Compare FSLE calculated from simulations with 
48 min output and 20 min output.
'''

import netCDF4 as netCDF
import tracpy
import tracpy.calcs
import numpy as np

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)

## Do everything for the 20 minute output case ##

d20 = netCDF.Dataset('tracks/2005-06-17T00gc.nc')
xg = d20.variables['xg'][50000:60000,:] # just use a subset
yg = d20.variables['yg'][50000:60000,:]
tp = d20.variables['tp'][:]
ind = xg==-1
lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
lonp[ind] = np.nan
latp[ind] = np.nan
xp[ind] = np.nan
yp[ind] = np.nan
d20.close()

# Calculate pairs
# let the index in axis 0 be the drifter id
ID = np.arange(lonp.shape[0])

# save pairs to save time since they are always the same
if not os.path.exists('tracks/pairs_limited.npz'):

    dist = np.zeros((lonp.shape[0],lonp.shape[0]))
    for idrifter in xrange(lonp.shape[0]):
        # dist contains all of the distances from other drifters for each drifter
        dist[idrifter,:] = tracpy.calcs.get_dist(lonp[idrifter,0], lonp[:,0], latp[idrifter,0], latp[:,0])
    pairs = []
    for idrifter in xrange(lonp.shape[0]):
        ind = find(dist[idrifter,:]<=1)
        for i in ind:
            if ID[idrifter] != ID[i]:
                pairs.append([min(ID[idrifter], ID[i]), 
                                max(ID[idrifter], ID[i])])

    pairs_set = set(map(tuple,pairs))
    pairs = map(list,pairs_set)# now pairs has only unique pairs of drifters
    pairs.sort() #unnecessary but handy for checking work
    np.savez('tracks/pairs_limited.npz', pairs=pairs)
else:
    pairs = np.load('tracks/pairs_limited.npz')['pairs']

# Run through pairs to save overall FSLE
fname = 'fsle_20min.npz'

dSave = np.zeros(20)
tSave = np.zeros(20)
nnans = np.zeros(20) # to collect number of non-nans over all drifters for a time
for ipair in xrange(len(pairs)):

    dSavetemp, tSavetemp = tracpy.calcs.calc_fsle(lonp[pairs[ipair][0],:], latp[pairs[ipair][0],:], 
                                lonp[pairs[ipair][1],:], latp[pairs[ipair][1],:], tp)
    ind = ~np.isnan(tSavetemp)
    dSave[ind] += dSavetemp[ind]
    tSave[ind] += tSavetemp[ind]
    nnans[ind] += 1
    # fsle += fsletemp
    # nnans += nnanstemp

# Save fsle for each file/area combination, NOT averaged
np.savez(fname, dSave=dSave, tSave=tSave, nnans=nnans)
print 'saved file', fname

## End 20 minute output case ##


## Do everything for the 48 minute output case ##

d48 = netCDF.Dataset('tracks/2005-06-17T00gc-N5_nsteps25.nc')
xg = d48.variables['xg'][50000:60000,:] # just use a subset
yg = d48.variables['yg'][50000:60000,:]
tp = d48.variables['tp'][:]
ind = xg==-1
lonp, latp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
lonp[ind] = np.nan
latp[ind] = np.nan
xp[ind] = np.nan
yp[ind] = np.nan
d48.close()

# Run through pairs to save overall FSLE
fname = 'fsle_48min.npz'

dSave = np.zeros(20)
tSave = np.zeros(20)
nnans = np.zeros(20) # to collect number of non-nans over all drifters for a time
for ipair in xrange(len(pairs)):

    dSavetemp, tSavetemp = tracpy.calcs.calc_fsle(lonp[pairs[ipair][0],:], latp[pairs[ipair][0],:], 
                                lonp[pairs[ipair][1],:], latp[pairs[ipair][1],:], tp)
    ind = ~np.isnan(tSavetemp)
    dSave[ind] += dSavetemp[ind]
    tSave[ind] += tSavetemp[ind]
    nnans[ind] += 1
    # fsle += fsletemp
    # nnans += nnanstemp

# Save fsle for each file/area combination, NOT averaged
np.savez(fname, dSave=dSave, tSave=tSave, nnans=nnans)
print 'saved file', fname

## End 48 minute output case ##

## Plots from the analysis ##
alpha = np.sqrt(2)
Rs = np.asarray([np.array([0.7])*alpha**i for i in np.arange(20)])

f = np.load('fsle_20min.npz')
l20 = 1./(f['tSave']/f['nnans'])
f.close()
f = np.load('fsle_48min.npz')
l48 = 1./(f['tSave']/f['nnans'])
f.close()

figure()
loglog(Rs[:-3], l20[:-3], 'go')
loglog(Rs[:-3], l48[:-3], 'bo')
