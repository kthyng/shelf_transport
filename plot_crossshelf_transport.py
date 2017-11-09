'''
Plot drifters from winter and summer transport regions, interannual. For paper.
Plot aggregated.

Example usage:
python3 plot_crossshelf_transport.py "summer" 'calc' --whichcalc '2Dtransport' --year 2004
python3 plot_crossshelf_transport.py "summer" 'calc' --whichcalc '1Dcrossing' --year 2004
python3 plot_crossshelf_transport.py "winter" 'calc' --whichcalc '1Dcrossing' --year 2004
python3 plot_crossshelf_transport.py "summer" 'plot' --whichcalc '1Dcrossing'
'''

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import netCDF4 as netCDF
import pdb
import matplotlib as mpl
import cmocean.cm as cmo
import cartopy
ccrs = cartopy.crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import xarray as xr
from pyproj import Proj
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm
import os
import argparse
import shapely.geometry


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



figextent = [-98, -87.5, 22.8, 30.5]
figsize = (7, 7)
top, right, left, bottom =.96, .98, .15, .01
caxpos = [0.19, 0.79, 0.24, 0.02]  # colorbar axis position
datex, datey = 0.01, 0.82  # location of date on figure
datax, datay = 0.41, 0.97  # location of data note on figure

merc = ccrs.Mercator(central_longitude=-85.0)
pc = ccrs.PlateCarree()
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
hlevs = [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450]  # isobath contour depths


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
grid = xr.open_dataset(grid_filename)


def calc(year, whichcalc):
    inputs = {'proj': 'lcc', 'ellps': 'clrk66', 'datum': 'NAD27',
              'lat_1': 22.5, 'lat_2': 31.0, 'lat_0': 30, 'lon_0': -94,
              'x_0': 0, 'y_0': 0}
    proj = Proj(**inputs)

    dataextent = [grid.lon_u.data.min(), grid.lon_u.data.max(), grid.lat_v.data.min(), grid.lat_v.data.max()]  # for histogram

    # set up for converting from xg/yg to lonp/latp:
    imt = grid.h.shape[1]  # 191
    jmt = grid.h.shape[0]  # 671
    X, Y = np.meshgrid(np.arange(imt), np.arange(jmt))
    # Triangulation for grid space to curvilinear space
    pts = np.column_stack((X.flatten(), Y.flatten()))
    tess = Delaunay(pts)
    tri = mtri.Triangulation(X.flatten(), Y.flatten(), tess.simplices.copy())
    flon = mtri.LinearTriInterpolator(tri, grid.lon_rho.data.flatten())
    flat = mtri.LinearTriInterpolator(tri, grid.lat_rho.data.flatten())
    fhg = mtri.LinearTriInterpolator(tri, grid.h.data.flatten())

    # set up projection for correct distances
    aecart = cartopy.crs.AzimuthalEquidistant(central_longitude=-96, central_latitude=28)
    ae = Proj(aecart.proj4_init)  # initialize proj4 projection for distances

    # how to choose drifters from which region to plot
    d = np.load('calcs/xyg0.npz')
    xg0 = d['xg0']; yg0 = d['yg0']
    lonp0 = flon(xg0+0.5, yg0+0.5)
    latp0 = flat(xg0+0.5, yg0+0.5)
    fh = mtri.LinearTriInterpolator(tri, grid.h.data.flatten())
    depths = fh(xg0, yg0)
    d.close()

    ishallow = depths < shelf_depth
    ideep = depths > shelf_depth

    # spacing for drifters
    if shelf_depth == 100:
        ds = 1  # 500
        ishelf_depth = 2
    elif shelf_depth == 50:
        ds = 250
        ishelf_depth = 1
    elif shelf_depth == 20:
        ds = 100
        ishelf_depth = 0

    # Set up isobaths for comparison
    if whichcalc == '1Dcrossing':
        # contour in lon/lat coords
        # cs = plt.contour(grid.lon_rho, grid.lat_rho , grid.h, [100])
        # contour in grid space
        cs = plt.contour(X, Y, grid.h, [100])
        p = cs.collections[0].get_paths()[0]  # pull out main isobath (not islands)
        v = p.vertices
        Iso = shapely.geometry.LineString(zip(v[:,0], v[:,1]))  # isobath as Line
        # save as lat lon too
        Iso_lon = flon(v[:,0], v[:,1])
        Iso_lat = flat(v[:,0], v[:,1])
        # save ll as Line
        Isoll = shapely.geometry.LineString(zip(Iso_lon, Iso_lat))  # isobath as Line
        # convert to equal distance projection
        Isoxy = aecart.project_geometry(Isoll, pc)[0]
        # Isopts = aecart.transform_points(pc, Iso_lon, Iso_lat)
        # Iso_x = Isopts[:,0]; Iso_y = Isopts[:,1]




    # Calculate the transport for each year for the winter wind transport area.
    # Want places where the depth is less than 100 meters and west of 90 deg.
    if region == 1: # shelf bend
        iwind = ideep * (lonp0<-95) * (lonp0>-97) * (latp0>26) * (latp0<28)
    elif region == 2: # Mississippi region
        iwind = ishallow * (lonp0<-88) * (lonp0>-89.3) * (latp0>28.7) * (latp0<30)
    elif region == 3: # winter wind transport region
        iwind = ishallow * (lonp0<-90) * (latp0>27.5)
    elif region == 4: # summer river transport region
        iwind = ishallow * (lonp0<-90) * (latp0>27.5)

    if whichcalc == '2Dtransport':
        fname = 'calcs/shelfconn/crossshelftransport-' + whichseason + '-' + str(year) + '.npz'
    elif whichcalc == '1Dcrossing':
        fname = 'calcs/shelfconn/crossshelfcrossing-' + whichseason + '-' + str(year) + '.npz'


    if not os.path.exists(fname):
        # initialize
        H = np.zeros((190, 670))
        # loop through tracks and add trajectories from region to plot
        if whichseason == 'winter':
            Files = glob('tracks/' + str(year) + '-0[1-2]-??T0?gc.nc')
        elif whichseason == 'summer':
            Files = glob('tracks/' + str(year) + '-0[7-8]-??T0?gc.nc')

        distances = []
        for File in Files:
            print(File)
            d = netCDF.Dataset(File)
            xg = d.variables['xg'][:]; yg = d.variables['yg'][:]
            tp = d['tp']

            # convert from xg/yg to xp/yp
            # nanind = (xg[iwind,:][::ds,:]==-1)
            # # Need to shift indices to move to rho grid of interpolator from
            # # arakawa c grid
            # xp = fx(xg[iwind,:][::ds,:]+0.5, yg[iwind,:][::ds,:]+0.5)
            # yp = fy(xg[iwind,:][::ds,:]+0.5, yg[iwind,:][::ds,:]+0.5)
            # xp[nanind] = np.nan; yp[nanind] = np.nan

            # convert from xg/yg to lonp/latp
            # nanind = (xg[iwind,:][::ds,:]==-1)
            # # Need to shift indices to move to rho grid of interpolator from
            # # arakawa c grid
            # lonp = flon(xg[iwind,:][::ds,:]+0.5, yg[iwind,:][::ds,:]+0.5)
            # latp = flat(xg[iwind,:][::ds,:]+0.5, yg[iwind,:][::ds,:]+0.5)
            # lonp[nanind] = np.nan; latp[nanind] = np.nan

            # indices of those that do or do not cross shelf, as desired
            shelffile = 'calcs/shelfconn/' + File.split('/')[-1][:-3] + '.npz'  # file of cross-shelf calculations
            # indices of drifters that cross 100 m isobath
            if cross:
                icross = ~np.isnan(np.load(shelffile)['cross'][ishelf_depth, iwind][::ds])
            else:
                icross = np.isnan(np.load(shelffile)['cross'][ishelf_depth, iwind][::ds])


            # rename xg/yg limited with indices by region and by crossings
            xgt = xg[iwind,:][icross,:]; ygt = yg[iwind,:][icross,:]


            if whichcalc == '2Dtransport':
                # Make histogram of drifter locations on numerical grid
                # subtract 2 from jmt, imt because psi grid
                Htemp, xe, ye = np.histogram2d(xgt.flatten(), ygt.flatten(), bins=[imt-1,jmt-1], range=[[0, imt-2], [0, jmt-2]])
                H += Htemp.T
            elif whichcalc == '1Dcrossing':
                # import pdb; pdb.set_trace()
                lines = [shapely.geometry.LineString(zip(xgt[i,~np.isnan(xgt[i,:])], ygt[i,~np.isnan(ygt[i,:])])) for i in range(xgt.shape[0])]
                # Lines = shapely.geometry.MultiLineString(lines)
                # find intersection point(s) of dline with Iso, but just take first
                for line in lines:
                    # find intersection of drifter trajectory and isobath
                    inter = line.intersection(Iso)
                    # define the Point of intersection
                    if isinstance(inter, shapely.geometry.point.Point):
                        # convert to lat/lon
                        lon = flon(*inter.coords[0]).data
                        lat = flat(*inter.coords[0]).data
                        pt = shapely.geometry.Point(inter.coords[0])
                        # pts.append(inter.coords[0])
                    elif isinstance(inter, shapely.geometry.multipoint.MultiPoint):
                        lon = flon(*inter[0].coords[0]).data
                        lat = flat(*inter[0].coords[0]).data
                        # convert from lat lon to equal distance projection
                        x, y = aecart.transform_point(lon, lat, pc)
                        # create Point
                        pt = shapely.geometry.Point(x, y)
                        # pts.append(inter[0].coords[0])
                    # find the distance along Iso to the intersection Point pt
                    # with project
                    # distances stores along-isobath distance in projected space
                    distances.append(Isoxy.project(pt))
                # x, y = zip(*pts)
                    # import pdb; pdb.set_trace()

        lon_psi = grid.lon_psi; lat_psi = grid.lat_psi

        if whichcalc == '2Dtransport':
            np.savez(fname, H=H, lon_psi=grid.lon_psi, lat_psi=grid.lat_psi)
        elif whichcalc == '1Dcrossing':
            np.savez(fname, distances=distances, Iso=Iso, Isoll=Isoll, Isoxy=Isoxy)

    else:
        if whichcalc == '2Dtransport':
            d = np.load(fname)
            H = d['H']; lon_psi = d['lon_psi']; lat_psi = d['lat_psi']

    return H, lon_psi, lat_psi


def plot(whichcalc="2Dtransport"):

    years = np.arange(2004, 2015)

    if whichcalc == '2Dtransport':

        if cross:
            cmap = cmo.ice_r
        else:
            cmap = 'darkcyan'


        fig, axarr = plt.subplots(4,3)
        fig.set_size_inches(8.9, 11.5)
        fig.subplots_adjust(left=0.04, bottom=0.1, right=1.0, top=0.99, wspace=0.0001, hspace=0.05)
        for i in range(12):

            # Titles for subplots
            if i==11:
                ax.set_axis_off()
                continue
            else:
                ax = fig.add_subplot(4,3,i+1, projection=merc)
                gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
                # the following two make the labels look like lat/lon format
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                # gl.xlocator = mticker.FixedLocator([-105, -95, -85, -75, -65])  # control where the ticks are
                # gl.xlabel_style = {'size': 15, 'color': 'gray'}  # control how the tick labels look
                # gl.ylabel_style = {'color': 'red', 'weight': 'bold'}
                gl.xlabels_bottom = False  # turn off labels where you don't want them
                gl.ylabels_right = False
                ax.add_feature(land_10m, facecolor='0.8')
                ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
                ax.add_feature(states_provinces, edgecolor='0.2')
                ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='0.2')
                ax.set_extent(figextent, pc)

            H, lon, lat = calc(years[i], whichcalc)

            # plot
            # ind = np.isnan(H)
            # H = np.ma.masked_where(ind, H)
            # import pdb; pdb.set_trace()
            mappable = ax.pcolormesh(lon, lat, H, cmap=cmap, transform=pc, vmin=0, vmax=100000)#, norm=LogNorm(vmin=2, vmax=1000))
            plt.colorbar(mappable)
            # ax.plot(lonp[icross,:].T, latp[icross,:].T, color=color, alpha=0.2, lw=0.1, transform=pc)

        if cross:
            plt.savefig('figures/transport/' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-cross_lowres.png', bbox_inches='tight', dpi=100)
            plt.savefig('figures/transport/' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-cross.png', bbox_inches='tight', dpi=300)
        else:
            plt.savefig('figures/transport/' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-notcross_lowres.png', bbox_inches='tight', dpi=100)
            plt.savefig('figures/transport/' + whichseason + str(region) + '_depth' + str(shelf_depth) + '_drifters-notcross.png', bbox_inches='tight', dpi=300)

    elif whichcalc == '1Dcrossing':




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('whichseason', type=str, help='"summer" or "winter"')
    parser.add_argument('which', type=str, help='"calc" or "plot"')
    parser.add_argument('--whichcalc', type=str, help='"2Dtransport" or "1Dcrossing"')
    parser.add_argument('--year', type=int, help='What year to plot, only for calc')

    args = parser.parse_args()

    whichseason = args.whichseason
    which = args.which
    whichcalc = args.whichcalc
    year = args.year

    shelf_depth = 100
    # region = 3
    # whichseason = 'summer'
    cross = True  # to plot the ones that do or do not cross shelf

    if whichseason == 'winter':
        region = 3
    elif whichseason == 'summer':
        region = 4

    if which == 'calc':
        calc(year, whichcalc)  # using input year
    elif which == 'plot':
        plot()
