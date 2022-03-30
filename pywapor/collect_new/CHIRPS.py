import datetime
from pywapor.general.logger import log
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap

def default_vars():
    variables =  {
        "precip": [("time", "latitude", "longitude"), "p"],
            }
    return variables

def default_post_processors():
    post_processors = []
    return post_processors

def fn_func(product_name, tile):
    fn = f"{product_name}_temp.nc"
    return fn

def url_func(product_name, tile):
    url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/chirps20GlobalDailyP05.nc?"
    return url

def download(folder, latlim, lonlim, timelim, 
                variables = None, post_processors = None):
    product_name = 'CHIRPS'
    tiles = [None]
    coords = {"x": "longitude", "y": "latitude", "t": "time"}
    variables = default_vars()
    post_processors = default_post_processors()
    data_source_crs = get_crss("WGS84")
    parallel = False
    spatial_tiles = False
    un_pw = accounts.get("NASA")
    request_dims = False
    ds = opendap.download(folder, product_name, latlim, lonlim, timelim, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims)
    return ds

if __name__ == "__main__":

    import pywapor.collect.accounts as accounts

    folder = r"/Users/hmcoerver/Downloads"

    latlim = [36.9, 43.7]
    lonlim = [5.2, 17.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 1)]

    # CHIRPS.
    ds = download(folder, latlim, lonlim, timelim)

