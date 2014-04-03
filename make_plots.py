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
from glob import glob

mpl.rcParams['text.usetex'] = True
mpl.rcParams.update({'font.size': 26})



def calc_histograms():
	'''
	Calculate the histograms for the connectivity calculations.
	Calculate both the histogram of starting locations and of connected 
	trajectories.
	'''


# Which type of plot: 'weatherband', 'seasonal', 'interannual'
which = 'seasonal'

# Which calculations to plot
base = 'calcs/shelfconn/'
Files = glob(base + '*-07-*.npz')
Files.extend(base + '*-08-*.npz')
Files.extend(base + '*-01-*.npz')
Files.extend(base + '*-02-*.npz')

# Number of bins to use in histogram
bins = (60,60)

# Grid info
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)

# Loop through calculation files
for File in Files:

	# Read in connectivity info (previously calculated)
	d = np.load(File)

	# Set up plot
	fig = plt.figure(figsize=(11,10))
	ax = fig.add_subplot(111)
	tracpy.plotting.background(grid=grid, ax=ax)

	## Make cross-shelf connectivity plots ##
	# Code from ipython notebook in shelf_transport/connectivity
	# Except here need to actually make histograms -- so I can play around with plot parameters

	# Histogram of starting locations
	Hstart = np.histogram2d(d['xp0'], d['yp0'], bins=bins, 
				range=[[grid['xpsi'].min(), grid['xpsi'].max()], 
						[grid['ypsi'].min(), grid['ypsi'].max()]])

	# Histogram of crossing locations
	Hcross = np.histogram2d(d['xp0'], d['yp0'], bins=bins, 
				range=[[grid['xpsi'].min(), grid['xpsi'].max()], 
						[grid['ypsi'].min(), grid['ypsi'].max()]])



	# ax.plot(xp[::5,tind], yp[::5,tind], 'o', color='darkcyan')
	# plt.savefig('figures/2010-07-01T00/' + str(tind) + '.png', dpi=150)
	# plt.close(fig)

	d.close()