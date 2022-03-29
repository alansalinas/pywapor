"""
https://opendap.cr.usgs.gov/opendap/hyrax/SRTMGL1_NUMNC.003/contents.html
"""
import datetime
import rasterio
import os
import json
import pywapor
import rasterio.crs
import rasterio.warp
import rioxarray.merge
import xarray as xr
import numpy as np
from functools import partial
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from pywapor.general.logger import log
from shapely.geometry import shape
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap
from pywapor.general.processing_functions import save_ds, create_selection
from pywapor.general import bitmasks

folder = r"/Users/hmcoerver/Downloads/merra2"

latlim = [26.9, 33.7]
lonlim = [25.2, 37.2]
# latlim = [28.9, 29.7]
# lonlim = [30.2, 31.2]
timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 15)]

variables = None
post_processors = None

product_name = "M2I1NXASM.5.12.4"

def tiles_intersect(latlim, lonlim):
    with open(os.path.join(pywapor.collect_new.__path__[0], "SRTM30_tiles.geojson")) as f:
        features = json.load(f)["features"]
    aoi = Polygon.from_bounds(lonlim[0], latlim[0], lonlim[1], latlim[1])
    tiles = list()
    for feature in features:
        shp = shape(feature["geometry"])
        tile = feature["properties"]["dataFile"]
        if shp.intersects(aoi):
            n = int(tile[1:3])
            e = int(tile[4:7])
            tiles.append((n, e))
    return tiles