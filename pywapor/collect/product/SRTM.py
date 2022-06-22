"""
https://opendap.cr.usgs.gov/opendap/hyrax/SRTMGL1_NUMNC.003/contents.html
"""
import datetime
import os
import json
import pywapor.collect
from pywapor.general.processing_functions import open_ds
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from shapely.geometry import shape
import pywapor.collect.protocol.opendap as opendap

def tiles_intersect(latlim, lonlim):
    with open(os.path.join(pywapor.collect.__path__[0], "product/SRTM30_tiles.geojson")) as f:
        features = json.load(f)["features"]
    aoi = Polygon.from_bounds(lonlim[0], latlim[0], lonlim[1], latlim[1])
    tiles = list()
    for feature in features:
        shp = shape(feature["geometry"])
        tile = feature["properties"]["dataFile"]
        if shp.intersects(aoi):
            tiles.append(tile.split(".")[0])
    return tiles

def default_vars(product_name, req_vars = ["z"]):

    variables = {
        "30M": {
            "SRTMGL1_DEM": [("time", "lat", "lon"), "z"],
            "crs": [(), "spatial_ref"],
            }
    }

    req_dl_vars = {
        "30M": {
            "z": ["SRTMGL1_DEM", "crs"],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def drop_time(ds):
    return ds.isel(time=0).drop("time")

def default_post_processors(product_name, req_vars = ["z"]):
    
    post_processors = {
        "30M": {
            "z": [],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def fn_func(product_name, tile):
    fn = f"{product_name}_{tile}.nc"
    return fn

def url_func(product_name, tile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/SRTMGL1_NC.003/{tile}.SRTMGL1_NC.ncml.nc4?"
    return url

def download(folder, latlim, lonlim, product_name = "30M", req_vars = ["z"], variables = None, post_processors = None, **kwargs):
    folder = os.path.join(folder, "SRTM")

    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    spatial_buffer = True
    if spatial_buffer:
        dx = dy = 0.0002777777777777768
        latlim = [latlim[0] - dy, latlim[1] + dy]
        lonlim = [lonlim[0] - dx, lonlim[1] + dx]

    timelim = [datetime.date(2000, 2, 10), datetime.date(2000, 2, 12)]
    tiles = tiles_intersect(latlim, lonlim)
        
    coords = {"x": ["lon", lonlim], "y": ["lat", latlim], "t": ["time", timelim]}
    
    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    data_source_crs = None
    parallel = False
    spatial_tiles = True
    un_pw = accounts.get("NASA")
    request_dims = True

    ds = opendap.download(folder, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims)

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 11)]

    # fn = os.path.join(os.path.join(folder, "SRTM"), "30M.nc")
    # if os.path.isfile(fn):
    #     os.remove(fn)

    # SRTM.
    ds = download(folder, latlim, lonlim)
    print(ds.rio.crs, ds.rio.grid_mapping)

