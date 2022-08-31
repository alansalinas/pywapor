import os
import pywapor.collect.protocol.sentinelapi as sentinelapi
from pywapor.general.curvilinear import regrid, create_grid
from pywapor.general.logger import log
import glob
import xarray as xr
import numpy as np
from datetime import datetime as dt
from pywapor.general.processing_functions import open_ds

def default_vars(product_name, req_vars):

    variables = {
        "SL_2_LST___": {
                    "LST_in.nc": [(), "lst"],
                    "geodetic_in.nc": [(), "coords"],
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

    post_processors = {
        "SL_2_LST___": {
            "lst": []
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def time_func(fn):
    start_dtime = np.datetime64(dt.strptime(fn.split("_")[7], "%Y%m%dT%H%M%S"))
    end_dtime = np.datetime64(dt.strptime(fn.split("_")[8], "%Y%m%dT%H%M%S"))
    dtime = start_dtime + (end_dtime - start_dtime)/2
    return dtime

def process_s3(scene_folder, variables, bb = None, **kwargs):

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
    
    product_folder = os.path.join(folder, "SENTINEL3")

    fn = os.path.join(product_folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        ds = open_ds(fn)
        if np.all([x in ds.data_vars for x in req_vars]):
            return ds
        else:
            ds = ds.close()

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

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

    ds = sentinelapi.process_sentinel(scenes, variables, process_s3, 
                                        time_func, f"{product_name}.nc", post_processors, bb = bb)

    return ds

if __name__ == "__main__":

    ...
    
    # folder = r"/Users/hmcoerver/On My Mac/create_table"
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    # # timelim = ["2021-07-01", "2021-07-11"]
    # timelim = ["2022-07-01", "2022-07-03"]

    # product_name = 'SL_2_LST___'

    # req_vars = ["lst"]
    # post_processors = None
    # variables = None
    # extra_search_kwargs = {}

    # ds = download(folder, latlim, lonlim, timelim, product_name, 
    #             req_vars, variables = variables,  post_processors = post_processors)


