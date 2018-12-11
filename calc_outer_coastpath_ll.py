'''
Calculate outerpath of coastal boxes in lat/lon using basemap and python 2
Run on hafen in ipython2 alias
'''

from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path
from pyproj import Proj

inputs = {'projection': 'lcc', 'llcrnrlon': -98.5, 'llcrnrlat': 22.5,
          'urcrnrlon': -87.5, 'urcrnrlat': 31.0, 'lat_0': 30,
          'lon_0': -94, 'resolution': 'i', 'area_thresh': 0.}

proj = Basemap(**inputs)

# read in outerpath in projected coordinates
outerpathxy = np.load('calcs/coastpaths.npz')['outerpathxy'].item()

# convert to lat/lon
verts = proj(*outerpathxy.vertices.T, inverse=True)

# make path
outerpathll = Path(np.vstack(verts).T)

# save
np.savez('calcs/coastpath_outerll.npz', outerpathll=outerpathll)
