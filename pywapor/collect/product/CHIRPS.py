from pywapor.general.logger import log
from pywapor.collect.protocol.projections import get_crss
import pywapor.collect.protocol.opendap as opendap
import pywapor.collect.accounts as accounts
import os
import numpy as np
from pywapor.general.processing_functions import open_ds

def default_vars(product_name, req_vars = ["p"]):
    
    variables =  {
        "P05": {
            "precip": [("time", "latitude", "longitude"), "p"],
        }
    }

    req_dl_vars = {
        "P05": {
            "p": ["precip"],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars = ["p"]):

    post_processors = {
        "P05": {
            "p": [],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def fn_func(product_name, tile):
    fn = f"{product_name}_temp.nc"
    return fn

def url_func(product_name, tile):
    url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/chirps20GlobalDailyP05.nc?"
    return url

def download(folder, latlim, lonlim, timelim, product_name = "P05", req_vars = ["p"],
                variables = None, post_processors = None):
    folder = os.path.join(folder, "CHIRPS")

    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    tiles = [None]
    coords = {"x": ["longitude", lonlim], "y": ["latitude", latlim], "t": ["time", timelim]}
    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    timedelta = np.timedelta64(30, "m")
    data_source_crs = get_crss("WGS84")
    parallel = False
    spatial_tiles = False
    un_pw = accounts.get("NASA")
    request_dims = False
    ds = opendap.download(folder, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims,
                timedelta = timedelta)
    return ds

if __name__ == "__main__":

    import datetime

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 11)]

    product_name = "P05"
    req_vars = ["p"]
    # CHIRPS.
    ds = download(folder, latlim, lonlim, timelim)

