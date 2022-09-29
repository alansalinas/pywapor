from pywapor.general.logger import log
from pywapor.collect.protocol.projections import get_crss
import pywapor.collect.protocol.opendap as opendap
import pywapor.collect.accounts as accounts
import os
import xarray as xr
import numpy as np
from pywapor.general.processing_functions import open_ds, remove_ds, save_ds

def default_vars(product_name, req_vars = ["p"]):
    """Given a `product_name` and a list of requested variables, returns a dictionary
    with metadata on which exact layers need to be requested from the server, how they should
    be renamed, and how their dimensions are defined.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list, optional
        List of variables to be collected, by default ["p"].

    Returns
    -------
    dict
        Metadata on which exact layers need to be requested from the server.
    """
    
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
    """Given a `product_name` and a list of requested variables, returns a dictionary with a 
    list of functions per variable that should be applied after having collected the data
    from a server.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list, optional
        List of variables to be collected, by default ["p"].

    Returns
    -------
    dict
        Functions per variable that should be applied to the variable.
    """

    post_processors = {
        "P05": {
            "p": [],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def fn_func(product_name, tile):
    """Returns a filename for the `product_name`.

    Parameters
    ----------
    product_name : str
        Name of the product.
    tile : None
        Not used.

    Returns
    -------
    str
        Filename.
    """
    fn = f"{product_name}_temp.nc"
    return fn

def url_func(product_name, tile):
    """Returns a url at which to collect CHIRPS data.

    Parameters
    ----------
    product_name : None
        Not used.
    tile : None
        Not used.

    Returns
    -------
    str
        The url.
    """
    url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/chirps20GlobalDailyP05.nc?"
    return url

def download(folder, latlim, lonlim, timelim, product_name = "P05", req_vars = ["p"],
                variables = None, post_processors = None):
    """Download CHIRPS data and store it in a single netCDF file.

    Parameters
    ----------
    folder : str
        Path to folder in which to store results.
    latlim : list
        Latitude limits of area of interest.
    lonlim : list
        Longitude limits of area of interest.
    timelim : list
        Period for which to prepare data.
    product_name : str, optional
        Name of the product to download, by default "P05".
    req_vars : list, optional
        Which variables to download for the selected product, by default ["p"].
    variables : dict, optional
        Metadata on which exact layers need to be requested from the server, by default None.
    post_processors : dict, optional
        Functions per variable that should be applied to the variable, by default None.

    Returns
    -------
    xr.Dataset
        Downloaded data.
    """
    folder = os.path.join(folder, "CHIRPS")

    appending = False
    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        os.rename(fn, fn.replace(".nc", "_to_be_appended.nc"))
        existing_ds = open_ds(fn.replace(".nc", "_to_be_appended.nc"))
        if np.all([x in existing_ds.data_vars for x in req_vars]):
            existing_ds = existing_ds.close()
            os.rename(fn.replace(".nc", "_to_be_appended.nc"), fn)
            existing_ds = open_ds(fn)
            return existing_ds[req_vars]
        else:
            appending = True
            fn = os.path.join(folder, f"{product_name}_appendix.nc")
            req_vars = [x for x in req_vars if x not in existing_ds.data_vars]

    tiles = [None]
    coords = {"x": ["longitude", lonlim], "y": ["latitude", latlim], "t": ["time", timelim]}
    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items() if k in req_vars}

    timedelta = np.timedelta64(30, "m")
    data_source_crs = get_crss("WGS84")
    parallel = False
    spatial_tiles = False
    un_pw = accounts.get("NASA")
    request_dims = False
    ds_new = opendap.download(fn, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims,
                timedelta = timedelta)


    if appending:
        ds = xr.merge([ds_new, existing_ds])
        lbl = f"Appending new variables (`{'`, `'.join(req_vars)}`) to existing file."
        ds = save_ds(ds, os.path.join(folder, f"{product_name}.nc"), encoding = "initiate", label = lbl)
        remove_ds(ds_new)
        remove_ds(existing_ds)
    else:
        ds = ds_new

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

