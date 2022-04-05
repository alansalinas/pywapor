import datetime
import pandas as pd
import pywapor.collect.accounts as accounts
from pywapor.general.logger import log
from pywapor.collect_new.protocol.projections import get_crss
import pywapor.collect_new.protocol.opendap as opendap
import fnmatch
from pywapor.collect_new.protocol.requests import find_paths

def default_vars(product_name):
    vars = {
        "M2I1NXASM.5.12.4": {
                    "T2M": [("time", "lat", "lon"), "t_air_i"],
                    "U2M": [("time", "lat", "lon"), "u2m_i"],
                    "V2M": [("time", "lat", "lon"), "v2m_i"],
                    "QV2M": [("time", "lat", "lon"), "qv_i"],
                    "TQV": [("time", "lat", "lon"), "wv_i"],
                    "PS": [("time", "lat", "lon"), "p_air_i"],
                    "SLP": [("time", "lat", "lon"), "p_air_0_i"],
                        },
        "M2T1NXRAD.5.12.4": {
                    "SWGNT": [("time", "lat", "lon"), "ra_i"],
                        }
    }
    return vars[product_name]

def default_post_processors(product_name):
    post_processors = {
        "M2I1NXASM.5.12.4": [],
        "M2T1NXRAD.5.12.4": [],
    }
    return post_processors[product_name]

def fn_func(product_name, tile):
    fn = f"{product_name}_{tile.strftime('%Y%m%d')}.nc"
    return fn

def url_func(product_name, tile):

    def _filter(tag):
        tag_value = tag["href"]
        if tag_value[-5:] == ".html":
            tag_value = tag_value[:-5] 
        return tag_value
    # Find the existing tiles for the given year and month, this is necessary
    # because the version number (`\d{3}`) in the filenames is irregular.
    regex = r"MERRA2_\d{3}\..*\.\d{8}.nc4.html"
    url = f"https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/{product_name}/{tile.year}/{tile.month:02d}/contents.html"
    tile_names = find_paths(url, regex, filter = _filter)

    # Find which of the existing tiles matches with the date.
    fn_pattern = f"MERRA2_*.*.{tile.strftime('%Y%m%d')}.nc4"
    fn = fnmatch.filter(tile_names, fn_pattern)[0]

    # Create the final url.
    url = f"https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/{product_name}/{tile.year}/{tile.month:02d}/{fn}.nc4?"  
    return url

def download(folder, latlim, lonlim, timelim, product_name,
                variables = None, post_processors = None):
    tiles = pd.date_range(timelim[0], timelim[1], freq="D")
    coords = {"x": ["lon", lonlim], "y":["lat", latlim], "t": ["time", timelim]}
    if isinstance(variables, type(None)):
        variables = default_vars(product_name)
    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name)
    data_source_crs = get_crss("WGS84")
    parallel = True
    spatial_tiles = False
    un_pw = accounts.get("NASA")
    request_dims = True
    ds = opendap.download(folder, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims)
    return ds


if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/merra2"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 15)]

    # MERRA2.
    # product_name = "M2I1NXASM.5.12.4"
    product_name = "M2T1NXRAD.5.12.4"

    ds = download(folder, latlim, lonlim, timelim, product_name)

