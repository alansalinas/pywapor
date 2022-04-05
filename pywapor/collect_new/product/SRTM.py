"""
https://opendap.cr.usgs.gov/opendap/hyrax/SRTMGL1_NUMNC.003/contents.html
"""
import datetime
import os
import json
import pywapor
import rasterio.crs
import rasterio.warp
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from pywapor.general.logger import log
from shapely.geometry import shape
import pywapor.collect_new.opendap as opendap

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

def default_vars():

    vars = {
        "SRTMGL1_DEM": [("time", "lat", "lon"), "z"],
        "crs": [(), "spatial_ref"],
            }

    return vars

def default_post_processors():
    post_processors = []
    return post_processors

def fn_func(product_name, tile):
    fn = f"{product_name}_N{tile[0]:02d}E{tile[1]:03d}.nc"
    return fn

def url_func(product_name, tile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/SRTMGL1_NC.003/N{tile[0]:02d}E{tile[1]:03d}.SRTMGL1_NC.ncml.nc4?"
    return url

def download(folder, latlim, lonlim, variables = None, post_processors = None):
    product_name = "SRTM"
    timelim = [datetime.date(2000, 2, 10), datetime.date(2000, 2, 12)]
    tiles = tiles_intersect(latlim, lonlim)
    coords = {"x": ["lon", lonlim], "y": ["lat", latlim], "t": ["time", timelim]}
    variables = default_vars()
    post_processors = default_post_processors()
    data_source_crs = None
    parallel = True
    spatial_tiles = True
    un_pw = accounts.get("NASA")
    request_dims = True

    ds = opendap.download(folder, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims)

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/dl_test_new"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]

    # SRTM.
    ds = download(folder, latlim, lonlim)
