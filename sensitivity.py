'''
Find coastal boxes containing environmental sensitivity index
of 10.
'''

import geopandas as gpd
import numpy as np
import shapely.geometry as geometry
import matplotlib.pyplot as plt


# proj = tracpy.tools.make_proj('nwgom')
# grid = tracpy.inout.readgrid('../../grid.nc', proj)

loc = '/Volumes/COWLITZ/shelf_transport/'

# Read in ESI data
# tx = gpd.read_file('ESI/ESI.shp')
tx = gpd.read_file('/Users/kthyng/Documents/research/papers/coastal-impacts-MPB/revisions/figures/ESI/ESI.shp')
txll = tx.to_crs({'proj': 'latlong'})

# bool for if extra sensitive to oil spill
# http://response.restoration.noaa.gov/maps-and-spatial-data/example-shoreline-ranking.html
tx10 = ['10' in tx.ESI[i] for i in range(len(tx))]
# # change relevant shapes to shapely shapes
# shapestx = [geometry.Polygon(geo) for geo in txll[tx10].geometry]
# # louisiana
la = gpd.read_file('/Users/kthyng/Documents/research/papers/coastal-impacts-MPB/revisions/figures/Louisiana_2003_Shapefiles/AVPROJ/SHAPE/esip.shp')
la10 = ['10' in la.ESI[i] for i in range(len(la))]

# Read in coastal boxes
d = np.load('calcs/coastpaths.npz')
paths_orig = d['paths']
d.close()

# give coastal boxes a buffer to overlap with sensitivity regions
# change to shapely polygons
paths = [geometry.Polygon(path.vertices) for path in paths_orig]
# buffered paths
pathsb = [path.buffer(-0.002) for path in paths]

# search for buffered coast boxes that overlap with vulnerable shapes
# THIS NEXT

# pathsb.intersects(txll[tx10])
# for each buffered coast box, find all intersections
overlaps = [txll[tx10].intersects(pathb) for pathb in pathsb]

# plot if ESI region contains '10'
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
# for path, pathb in zip(paths, pathsb):
#     # show buffered box
#     ax.fill(*pathb.boundary.xy, color='0.4', alpha=1.0)
#     # show original coast box
#     ax.fill(*path.boundary.xy, color='0.1', alpha=1.0)
# # color coast box if overlaps with vulnerable shape
# for pathb, overlap in zip(pathsb, overlaps):
#     if overlap.sum() > 0:
#         ax.fill(*pathb.boundary.xy, color='r', alpha=0.7)
#         ax.fill(*path.boundary.xy, color='r', alpha=1.0)
# overlay ESI shapes
# tx.to_crs({'proj': 'latlong'}).plot(ax=ax, color='k')
txll.plot(ax=ax, color='k')
# color ESI shapes if vulnerable
txll[tx10].plot(ax=ax, color='blue', linewidth=2)

la.plot(ax=ax, color='k')
# color ESI shapes if vulnerable
la[la10].plot(ax=ax, color='blue', linewidth=2)
plt.show()

# 
