import os
import pywapor.collect.protocol.sentinelapi as sentinelapi
from pywapor.general.curvilinear import regrid, create_grid
from pywapor.general.logger import log, adjust_logger
import glob
import xarray as xr
import numpy as np
from datetime import datetime as dt
from pywapor.general.processing_functions import open_ds, remove_ds, save_ds

def default_vars(product_name, req_vars):
    """Given a `product_name` and a list of requested variables, returns a dictionary
    with metadata on which exact layers need to be requested from the server, how they should
    be renamed, and how their dimensions are defined.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list
        List of variables to be collected.

    Returns
    -------
    dict
        Metadata on which exact layers need to be requested from the server.
    """
    variables = {
        "SL_2_LST___": {
                    "LST_in.nc": [(), "lst", []],
                    "geodetic_in.nc": [(), "coords", []],
        }
    }

    req_dl_vars = {
        "SL_2_LST___": {
            "lst": ["LST_in.nc", "geodetic_in.nc"]
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out

def default_post_processors(product_name, req_vars):
    """Given a `product_name` and a list of requested variables, returns a dictionary with a 
    list of functions per variable that should be applied after having collected the data
    from a server.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list
        List of variables to be collected.

    Returns
    -------
    dict
        Functions per variable that should be applied to the variable.
    """
    post_processors = {
        "SL_2_LST___": {
            "lst": []
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def time_func(fn):
    """Return a np.datetime64 given a filename.

    Parameters
    ----------
    fn : str
        Filename.

    Returns
    -------
    np.datetime64
        Date as described in the filename.
    """
    start_dtime = np.datetime64(dt.strptime(fn.split("_")[7], "%Y%m%dT%H%M%S"))
    end_dtime = np.datetime64(dt.strptime(fn.split("_")[8], "%Y%m%dT%H%M%S"))
    dtime = start_dtime + (end_dtime - start_dtime)/2
    return dtime

def s3_processor(scene_folder, variables, bb = None, **kwargs):

    ncs = [glob.glob(os.path.join(scene_folder, "**", "*" + k), recursive = True)[0] for k in variables.keys()]

    ds = xr.open_mfdataset(ncs)

    ds = ds.set_coords(("longitude_in", "latitude_in"))
    ds = ds.rename_vars({"longitude_in": "x", "latitude_in": "y"})
    ds = ds.rename_dims({"rows": "ny", "columns": "nx"})
    ds = ds[["LST", "LST_uncertainty"]]

    ds = ds.where(ds.LST_uncertainty < 2.5)
    ds = ds.drop_vars("LST_uncertainty")

    grid_ds = create_grid(ds, 0.01, 0.01, bb = bb)
    ds = regrid(grid_ds, ds)
    ds = ds.rio.write_crs(4326)

    ds = ds.rename_vars({"LST": "lst"})

    return ds

def download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = None,  post_processors = None,
                extra_search_kwargs = {}):
    """Download SENTINEL3 data and store it in a single netCDF file.

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
    product_name : str
        Name of the product to download.
    req_vars : list
        Which variables to download for the selected product.
    variables : dict, optional
        Metadata on which exact layers need to be requested from the server, by default None.
    post_processors : dict, optional
        Functions per variable that should be applied to the variable, by default None.
    extra_search_kwargs : dict
        Extra search kwargs passed to SentinelAPI, by default {}.

    Returns
    -------
    xr.Dataset
        Downloaded data.
    """
    product_folder = os.path.join(folder, "SENTINEL3")

    appending = False
    fn = os.path.join(product_folder, f"{product_name}.nc")
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
            fn = os.path.join(product_folder, f"{product_name}_appendix.nc")
            req_vars = [x for x in req_vars if x not in existing_ds.data_vars]

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items() if k in req_vars}

    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    search_kwargs = {
                        "platformname": "Sentinel-3",
                        "producttype": product_name,
                        # "limit": 10,
                        }

    search_kwargs = {**search_kwargs, **extra_search_kwargs}

    # NOTE node_filter not supported yet for Sentinel-3, so downloading all files for now.
    # https://github.com/sentinelsat/sentinelsat/issues/566
    scenes = sentinelapi.download(product_folder, latlim, lonlim, timelim, 
                                    search_kwargs, node_filter = None)

    ds_new = sentinelapi.process_sentinel(scenes, variables, "SENTINEL3", time_func, os.path.split(fn)[-1], post_processors, bb = bb)

    if appending:
        ds = xr.merge([ds_new, existing_ds])
        lbl = f"Appending new variables (`{'`, `'.join(req_vars)}`) to existing file."
        ds = save_ds(ds, os.path.join(product_folder, f"{product_name}.nc"), encoding = "initiate", label = lbl)
        remove_ds(ds_new)
        remove_ds(existing_ds)
    else:
        ds = ds_new

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Local/s3_test"
    adjust_logger(True, folder, "INFO")
    timelim = ["2022-03-25", "2022-04-15"]
    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]

    product_name = 'SL_2_LST___'

    req_vars = ["lst"]
    post_processors = None
    variables = None
    extra_search_kwargs = {}

    ds = download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = variables,  post_processors = post_processors)


