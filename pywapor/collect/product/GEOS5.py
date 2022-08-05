from pywapor.collect.protocol import opendap
from pywapor.collect.protocol.projections import get_crss
import os
from pywapor.general.processing_functions import open_ds
from pywapor.enhancers.temperature import kelvin_to_celsius
from functools import partial
import numpy as np
from pywapor.enhancers.pressure import pa_to_kpa

def default_vars(product_name, req_vars):

    variables = {
        "inst3_2d_asm_Nx": {
                    "t2m": [("time", "lat", "lon"), "t_air"],
                    "u2m": [("time", "lat", "lon"), "u2m"],
                    "v2m": [("time", "lat", "lon"), "v2m"],
                    "qv2m": [("time", "lat", "lon"), "qv"],
                    "tqv": [("time", "lat", "lon"), "wv"],
                    "ps": [("time", "lat", "lon"), "p_air"],
                    "slp": [("time", "lat", "lon"), "p_air_0"],
                        },
    }

    req_dl_vars = {
        "inst3_2d_asm_Nx": {
            "t_air": ["t2m"],
            "t_air_max": ["t2m"],
            "t_air_min": ["t2m"],
            "u2m": ["u2m"],
            "v2m": ["v2m"],
            "qv": ["qv2m"],
            "wv": ["tqv"],
            "p_air": ["ps"],
            "p_air_0": ["slp"],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out

# def pa_to_kpa(ds, var):
#     ds[var] = ds[var] / 1000
#     return ds

def default_post_processors(product_name, req_vars):

    post_processors = {
        "inst3_2d_asm_Nx": {
            "t_air": [kelvin_to_celsius], 
            "t_air_max": [partial(kelvin_to_celsius, in_var = "t_air", out_var = "t_air_max")],
            "t_air_min": [partial(kelvin_to_celsius, in_var = "t_air", out_var = "t_air_min")],
            "u2m": [],
            "v2m": [],
            "qv": [],
            "wv": [],
            "p_air": [pa_to_kpa],
            "p_air_0": [pa_to_kpa],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def download(folder, latlim, lonlim, timelim, product_name, req_vars,
                 variables = None, post_processors = None):

    folder = os.path.join(folder, "GEOS5")
    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    spatial_buffer = True
    if spatial_buffer:
        latlim = [latlim[0] - 0.25, latlim[1] + 0.25]
        lonlim = [lonlim[0] - 0.3125, lonlim[1] + 0.3125]

    coords = {"x": ["lon", lonlim], "y": ["lat", latlim], "t": ["time", timelim]}

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    data_source_crs = get_crss("WGS84")

    url = f"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/{product_name}"

    fp = os.path.join(folder, f"{product_name}.nc")

    timedelta = np.timedelta64(90, "m")

    ds = opendap.download_xarray(url, fp, coords, variables, post_processors, 
                                    data_source_crs = data_source_crs,
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

    # GEOS5.
    product_name = "inst3_2d_asm_Nx"

    variables = None
    post_processors = None

    req_vars = [
                "t_air", "u2m", "v2m", "qv", 
                "wv", 
                "p_air", "p_air_0"
                ]

    ds = download(folder, latlim, lonlim, timelim, product_name, req_vars = req_vars)
    print(ds.rio.crs, ds.rio.grid_mapping)