'''
Taking code from plot_conn2coast.py to recreate the interannual summer plot
for DWH paper. Run this on tahoma.
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, cbook
import matplotlib.patches as Patches
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean.cm as cmo
import numpy as np
import xarray as xr
import os
import cartopy.mpl.geoaxes
import tabs

pc = cartopy.crs.PlateCarree()
merc = cartopy.crs.Mercator(central_longitude=-85)
# lcc is the cartopy version of the basemap projection I used to use
lcc = cartopy.crs.LambertConformal(central_longitude=-94.0, central_latitude=30,
                                   false_easting=466061.13903578784,
                                   false_northing=827157.5974829497)
land_10m = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cartopy.feature.COLORS['land'])
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
cmap = cmo.speed
levels = np.linspace(0,100,11)
grid = xr.open_dataset('../grid.nc')
season = 'summer' # 'winter' or 'summer'
whichboxes = 'galveston' # 'all', 'porta', 'galveston', 'both'
bins = (100,100)
levels = np.linspace(0,100,11)
extent = [-98, -87.5, 22.8, 30.5]
hlevs = [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450]  # isobath contour depths
# https://stackoverflow.com/questions/55385515/embed-small-map-cartopy-on-matplotlib-figure
axes_class = cartopy.mpl.geoaxes.GeoAxes
axes_kwargs=dict(map_projection=merc)

def resize(A, dim):
    """
    Average neighboring elements in an array A along a dimension dim.
    Args:
        A (array): array size, [m x n] in 2D case. Can be up to 3D.
        dim (int): dimension on which to act. Can be up to 2 (0, 1, 2).
    Returns:
        * A - array of size [(m-1) x n] if dim is 0
    """

    # B is A but rolled such that the dimension that is to be resized is in
    # the 0 position
    B = np.rollaxis(A, dim)

    # Do averaging
    B = 0.5*(B[0:-1]+B[1:])

    # Roll back to original
    return np.rollaxis(B, 0, dim+1)

# load in initial drifter starting locations in grid space
d = np.load('calcs/xyp0.npz', encoding='latin1')
xp0 = d['xp0']; yp0 = d['yp0']
d.close()

# # Calculate xrange and yrange for histograms
# Xrange = [grid.x_psi.min(), grid.x_psi.max()]
# Yrange = [grid.y_psi.min(), grid.y_psi.max()]

# Save a histogram of the number of drifters that started in each bin
Hstartfile = 'calcs/coastconn/likelihood/Hstart.npz'
if not os.path.exists(Hstartfile):
    Hstart, xe, ye = np.histogram2d(xp0, yp0, bins=bins,
                range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
    np.savez(Hstartfile, Hstart=Hstart, xe=xe, ye=ye)
else:
    d = np.load(Hstartfile, encoding='latin1')
    Hstart = d['Hstart']
    xe = d['xe']; ye = d['ye']
    d.close()


# Find indices of all drifters that start in the coastal boxes
# start indices of drifters that start in each box path
pts = np.load('calcs/alongcoastconn/inds-in-coast-paths.npz', encoding='latin1')['pts']

# Port A:
if whichboxes=='all':
    boxes = np.arange(0,342)
    whichH = 'Hall' # Hall for histogram for entire coastline at once; H for by coast box
    zoomed = False
elif whichboxes=='all2':
    boxes = np.arange(0,342)
    whichH = 'H' # use coast box histograms instead of combined
    zoomed = False
elif whichboxes=='pacyto':  # port aransas cytobot
    boxes = np.arange(112,115)
    whichH = 'H' # use coast box histograms instead of combined
    zoomed = False
    plume = False
elif whichboxes=='sscyto':  # surfside cytobot
    boxes = np.arange(155,158)
    whichH = 'H' # use coast box histograms instead of combined
    zoomed = False
    plume = False
elif whichboxes=='bothcyto':  # pa and surfside cytobot
    boxes = np.hstack((np.arange(112,115), np.arange(155,158)))
    whichH = 'H' # use coast box histograms instead of combined
    zoomed = False
    plume = False
elif whichboxes=='porta':
    boxes = np.arange(103,123)
    whichH = 'H'
    zoomed = True
    # limits for zoomed box
    if season == 'winter':
    # x1, x2, y1, y2 = 86000, 340800, 465000, 715000
        x1, x2, y1, y2 = 86000, 277100, 527500, 715000
    elif season == 'summer':
        x1, x2, y1, y2 = 86000, 277100, 465000, 652500
    zoom = 2.5
    plume = False
elif whichboxes=='galveston':
    boxes = np.arange(160,180)
    whichH = 'H'
    zoomed = True
    x1, x2, y1, y2 = 277100, 531900, 660000, 810000
    zoom = 2.75
    plume = False
    # x1, x2, y1, y2 = 277100, 531900, 560000, 810000
    # zoom = 2.0
elif whichboxes=='both':
    boxes = np.arange(103,180)
    whichH = 'H'
    zoomed = True
    x1, x2, y1, y2 = 86000, 531900, 465000, 810000
    zoom = 1.5

d = np.load('calcs/coastpaths.npz', encoding='latin1')
pathsxy = d['pathsxy']
d.close()

## Read in files ##
filename = 'calcs/coastconn/likelihood/hist-' + season + 'interannual-' + whichboxes + '.npz'
if not os.path.exists(filename):
    base = 'calcs/coastconn/likelihood/hist-'
    Files = []
    if season == 'winter':
        months = '[1-2]'
    elif season == 'summer':
        months = '[7-8]'
    Files.append(glob.glob(base + '2004-0' + months + '.npz'))
    Files.append(glob.glob(base + '2005-0' + months + '.npz'))
    Files.append(glob.glob(base + '2006-0' + months + '.npz'))
    Files.append(glob.glob(base + '2007-0' + months + '.npz'))
    Files.append(glob.glob(base + '2008-0' + months + '.npz'))
    Files.append(glob.glob(base + '2009-0' + months + '.npz'))
    Files.append(glob.glob(base + '2010-0' + months + '.npz'))
    Files.append(glob.glob(base + '2011-0' + months + '.npz'))
    Files.append(glob.glob(base + '2012-0' + months + '.npz'))
    Files.append(glob.glob(base + '2013-0' + months + '.npz'))
    Files.append(glob.glob(base + '2014-0' + months + '.npz'))
    H = np.zeros((len(Files),6,100,100))
    ndbox = np.zeros((len(Files),6,342))
    for i,files in enumerate(Files): # winter and summer
        numfiles = 0
        for File in files: # months/years within winter or summer
            print(File)
            d = np.load(File, encoding='latin1')
            numfiles += d['numfiles']
            idstemp = d['ids'].item()
            for k in xrange(d['numfiles']): # loop over simulations in file which don't want to mix up
                for j in xrange(6): #advection days
                    # add together the drifters for the desired boxes
                    ids = []
                    [ids.extend(idstemp[k,j,box]) for box in boxes]
                    ids = list(set(ids)) # eliminate nonunique drifters
                    Htemp, _, _ = np.histogram2d(xp0[ids], yp0[ids], bins=bins,
                                    range=[[Xrange[0], Xrange[1]], [Yrange[0], Yrange[1]]])
                    H[i,j] += Htemp
            days = d['days']
            xe = d['xe']; ye = d['ye']
            ndbox[i,:,:] += d['ndbox']
            d.close()
        # Divide by number of starting drifters
        H[i,:,:,:] /= (numfiles*Hstart)
    np.savez(filename, H=H, days=days, xe=xe, ye=ye, ndbox=ndbox)

else:
    d = np.load(filename, encoding='latin1')
    H = d['H']; days = d['days']; xe = d['xe']; ye = d['ye']; ndbox = d['ndbox']
    d.close()
####

# limit boxes
pathsxy = pathsxy[boxes]
ndbox = ndbox[:,:,boxes]

XE, YE = np.meshgrid(resize(xe, 0), resize(ye, 0))
# xe, ye = resize(xe, 0), resize(ye, 0)

# Loop through advection days
# for j in xrange(H.shape[1]):
j = H.shape[1]-1

## Plot setup ##
def setup(ax, left=False, top=False):
    ax.set_extent(extent, pc)
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_bottom = False  # turn off labels where you don't want them
    gl.ylabels_right = False
    gl.xlocator = ticker.FixedLocator(np.arange(-98, -86, 2))
    gl.ylocator = ticker.FixedLocator(np.arange(24, 32, 2))
    if not left:
        gl.ylabels_left = False
    if not top:
        gl.xlabels_top = False
    ax.add_feature(land_10m, facecolor='0.8')
    ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
    ax.add_feature(states_provinces, edgecolor='0.2')
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', edgecolor='0.2')
    # plot isobaths
    ax.contour(grid.lon_rho, grid.lat_rho, grid.h, hlevs, colors='0.6',
               transform=pc, linewidths=0.5, alpha=0.5)

fig, axarr = plt.subplots(4,3, subplot_kw=dict(projection=merc))
fig.set_size_inches(8.7, 10.5)
fig.subplots_adjust(left=0.008, bottom=0.08, right=1.0, top=0.98, wspace=0.07, hspace=0.005)
for i, ax in enumerate(axarr.flatten()):
    if i in [3,6,9]:
        setup(ax, left=True, top=False)
    elif i == 0:
        setup(ax, left=True, top=True)
    elif i in [1,2]:
        setup(ax, left=False, top=True)
    elif i==11:
        ax.outline_patch.set_visible(False)
        continue
    else:
        setup(ax)
    if i!=11:
        ax.outline_patch.set_visible(False)
        var = H[i,j].T*100.
        ax.text(0.07, 0.88, str(2004+i), transform=ax.transAxes, fontsize=14)
        mappable = ax.contourf(XE, YE, var, cmap=cmap, levels=levels, transform=lcc)


        # Add on vulnerability of coastline
        # Need to plot the coast boxes as patches and color them according to vulnerability level
        # http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
        # we need to normalize the data to 0..1 for the full
        # range of the colormap
        fracs = ndbox[i,j,:].astype(float)/ndbox[:,j,:].max() # max across years
        norm = colors.Normalize(fracs.min(), fracs.max())

        # Save patches together
        patches = []
        for path in pathsxy:
            patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5, transform=lcc))

        # assign shades of colormap to the patches according to values, and plot
        for thisfrac, thispatch in zip(fracs, patches):
            color = cmo.matter(norm(thisfrac))
            thispatch.set_facecolor(color)
            ax.add_patch(thispatch)

        if zoomed:# and j==H.shape[1]-1: # magnification for longest advection time available

            # Save patches together
            patches = []
            for path in pathsxy:
                patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5, transform=lcc))
                # ax.add_patch(patch)

            # Inset image
            axins = zoomed_inset_axes(ax, zoom, loc=4, axes_class=axes_class, axes_kwargs=axes_kwargs) # zoom=6
            setup(axins)
            axins.contourf(XE, YE, var, cmap=cmap, levels=levels, transform=lcc)

            # assign shades of colormap to the patches according to values, and plot
            for thisfrac, thispatch in zip(fracs, patches):
                color = cmo.matter(norm(thisfrac))
                thispatch.set_facecolor(color)
                axins.add_patch(thispatch)

            # subregion of the original image
            axins.set_extent([x1, x2, y1, y2], lcc)
            # axins.set_xlim(x1,x2)
            # axins.set_ylim(y1,y2)
            plt.xticks(visible=False)
            plt.yticks(visible=False)
            plt.setp(axins,xticks=[],yticks=[])
            # draw a bbox of the region of the inset axes in the parent axes and
            # connecting lines between the bbox and the inset axes area
            mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")#, kwargs=dict(zorder=8))
            plt.draw()
            plt.show()

            # overlay wind and salinity data
            buoy = '42035'
            df = tabs.read(buoy, str(2004+i) + '-7-15', str(2004+i) + '-9-15', model=True)
            Ukey = '%s: East [m/s] (air)' % buoy
            Vkey = '%s: North [m/s] (air)' % buoy

            # # check there are enough data points
            # assert df[Ukey].count()/len(df) > 0.9

            # plot mean wind
            axins.quiver(np.array([-95.7]), np.array([29.5]),
                         np.array([df[Ukey].mean()]), np.array([df[Vkey].mean()]),
                         transform=pc, color='k', alpha=0.6, scale=9, width=0.03,
                         pivot='middle', zorder=10, scale_units='width')
# , headlength=3, headaxislength=2.8
            # plot mean subtidal salt
            key = '%s: Salinity' % buoy

            # # check there are enough data points
            # assert df[key].count()/len(df) > 0.9
            # model has hourly output
            meansalt = df[key].rolling(center=True, window=24).mean().mean()
            axins.scatter([-95.3], [29.5], c=[meansalt], s=100, cmap=cmo.haline,
                          transform=pc, vmin=20, vmax=34, zorder=10, edgecolors='k')
            axins.text(-95.4, 29.25, '%2.0f' % meansalt, transform=pc)

            # ax.set_frame_on(False)
cax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) #colorbar axes
cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
cb.set_label('Likelihood of hitting the coast in ' + str(days[j]) + ' days [%]')
if zoomed:
    fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + 'zoomed_lowres.png', bbox_inches='tight')
    fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + 'zoomed.png', bbox_inches='tight', dpi=300)
else:
    fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + '_lowres.png', bbox_inches='tight')
    fig.savefig('figures/coastconn/likelihood/interannual-' + season + str(days[j]) + 'days' + whichboxes + '.png', bbox_inches='tight', dpi=300)
# plt.close()
