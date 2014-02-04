'''
Run plots for drifters from these simulations.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob

mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 26})

# Tracks to plot
File = 'tracks/2010-07-01T00gc.nc'

# Read in info
d = netCDF.Dataset(File)
xg = d.variables['xg'][:]
yg = d.variables['yg'][:]
tg = d.variables['tp'][:]
d.close()

# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)

# Change to projected drifter locations now
nanind = np.isnan(xg) # indices where nans are location in xg, yg; for reinstitution of nans
xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
xp[nanind] = np.nan; yp[nanind] = np.nan
del(xg,yg) # don't need grid info anymore

# Loop through time
for tind in xrange(2):#tg.size):

	# Set up plot
	fig = plt.figure(figsize=(11,10))
	ax = fig.add_subplot(111)
	tracpy.plotting.background(grid=grid, ax=ax)

	ax.plot(xp[::5,tind], yp[::5,tind], 'o', color='darkcyan')
	plt.savefig('figures/2010-07-01T00/' + str(tind) + '.png', dpi=150)
	plt.close(fig)
