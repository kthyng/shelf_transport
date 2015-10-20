'''
Plot interannual standard deviation in cross-shelf transport, for winter and summer.
'''


import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
import matplotlib as mpl
import numpy as np
import op


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


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

fig, axarr = plt.subplots(1,2)#, sharex=True)
fig.set_size_inches(13, 6.6125)
fig.subplots_adjust(left=0.045, bottom=0.15, right=1.0, top=0.96, wspace=0.005, hspace=0.04)

cmap = 'Purples'
for i, ax in enumerate(axarr):

   # Titles for subplots
    if i==0:
        tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
        ax.set_title('Winter')
        d = np.load('figures/cross/interannual-winter100H.npz')
        XE, YE = np.meshgrid(op.resize(d['xe'], 0), op.resize(d['ye'], 0))
        mappable = ax.contourf(XE, YE, np.std(d['H'], axis=0).T, cmap=cmap, vmin=0, vmax=32)#, levels=levels, extend=extend)
        ax.contour(grid['xr'], grid['yr'], grid['h'], [100], colors='0.1', linewidth=3)

    elif i==1:
        tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
        ax.set_title('Summer')
        d = np.load('figures/cross/interannual-summer100H.npz')
        mappable = ax.contourf(XE, YE, np.std(d['H'], axis=0).T, cmap=cmap, vmin=0, vmax=32)#, levels=levels, extend=extend)
        ax.contour(grid['xr'], grid['yr'], grid['h'], [100], colors='0.1', linewidth=3)

    ax.set_frame_on(False)

    # Colorbar
    cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
    cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')

    cb.set_label('Probability of drifters crossing shelf (%)')

fig.savefig('figures/cross/seasonal-std.png', bbox_inches='tight')
fig.savefig('figures/cross/seasonal-std-highres.png', bbox_inches='tight', dpi=300)