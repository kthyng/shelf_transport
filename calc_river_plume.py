'''
Looking at how to define river plume.
'''


import netCDF4 as netCDF
import os

d = netCDF.Dataset('tracks/2007-08-01T00gc.nc')

xg = d['xg'][:]
yg = d['yg'][:]
tp = d['tp'][:]

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
nc = netCDF.Dataset(loc)

# number of time indices to perform analysis over
days = (tp - tp[0])/86400
ndays = 30  # number of days for analysis
ntimes = np.where(days <= ndays)[0][-1] + 1

salt = tracpy.calcs.Var(xg, yg, tp, 'salt', nc)

# indices of drifters that have salinity always below psu
ilt = {}
for psu in np.arange(32,36):
    # make sure it is below for all times
    ilt[psu] = (salt[:,:ntimes]<=psu).sum(axis=1) == ntimes

# indices of drifters that have salinity always above psu
igt = {}
for psu in np.arange(32,36):
    igt[psu] = (salt[:,:ntimes]>psu).sum(axis=1) == ntimes

# indices of drifters that go from below to above or from above to below
iltgt = {}
igtlt = {}
for psu in np.arange(32,36):
    # make sure it is both below and above, and starts below
    iltgt[psu] = ((salt[:,:ntimes]<=psu).sum(axis=1) < ntimes) & (salt[:,0]<=psu)
    igtlt[psu] = ((salt[:,:ntimes]<=psu).sum(axis=1) < ntimes) & (salt[:,0]>psu)




os.makedirs('figures/salt_dist', exist_ok=True)

for psu in np.arange(32,36):
    # overall salinity
    figure(figsize=(8,4)); plt.hist(salt[:,:ntimes].flatten(), bins=36, range=(0,36), alpha=0.3)

    # always less than psu
    plt.hist(salt[ilt[psu],:ntimes].flatten(), bins=36, range=(0,36), alpha=0.3)

    # always greater than psu
    plt.hist(salt[igt[psu],:ntimes].flatten(), bins=36, range=(0,36), alpha=0.3)

    # less than, then greater than
    #plt.hist(salt[iltgt[psu],:ntimes].flatten(), bins=36, range=(0,36), alpha=0.3)

    # greater than, then less than
    #plt.hist(salt[igtlt[psu],:ntimes].flatten(), bins=36, range=(0,36), alpha=0.3)

    plt.title('psu=%i' % psu)
    plt.xlabel('Salinity [psu]')
    plt.ylabel('Number of drifters in bin')

    plt.savefig('figures/salt_dist/%i_ndays%i.png' % (psu,ndays), bbox_inches='tight')
    plt.savefig('figures/salt_dist/%i_ndays%i.pdf' % (psu,ndays), bbox_inches='tight')
    plt.close(plt.gcf())




## look at salinity by salinity class and time

# ndays = 5  # number of days for analysis
ndays = np.arange(1,11)
ds = 0.5
for nday in ndays:
    ntimes = np.where(days <= nday)[0][-1] + 1

    #psu = 30
    nums = []
    psus = np.arange(29,36,ds)
    for psu in psus:
        # indices of drifters that are in band of psu for all ntimes
        inds = ((salt[:,:ntimes] >= psu) & (salt[:,:ntimes] < psu+ds))
        nums.append(((inds.sum(axis=1) == ntimes).sum()/inds[:,0].sum(axis=0).sum())*100)

    plt.plot(psus, nums, label='%i days' % nday)

plt.xlabel('Salinity [psu]')
plt.ylabel('% drifters that stay in salinity class')
plt.legend(loc='best')
plt.savefig('figures/salt_dist/mixing.pdf', bbox_inches='tight')
plt.savefig('figures/salt_dist/mixing.png', bbox_inches='tight')
