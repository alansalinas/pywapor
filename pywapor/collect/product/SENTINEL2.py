import glob
import os
import xarray as xr
import numpy as np
from pywapor.collect.product.Landsat.C2L2SP import calc_albedo, calc_ndvi, mask_invalid, open_as_xr
from datetime import datetime as dt
from pywapor.general.logger import log
import pywapor.collect.protocol.sentinelapi as sentinelapi
import numpy as np
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.processing_functions import open_ds

def process_s2(scene_folder, variables):

    fps = {v[1]: glob.glob(os.path.join(scene_folder, "**", "*" + k), recursive = True)[0] for k, v in variables.items()}

    ds = xr.concat([open_as_xr(fp, name) for name, fp in fps.items()], "band")

    qa_da = ds.band_data.sel(band = "qa")
    data = ds.where(ds.band != "qa", drop = True)

    valid_range = {
            "blue":     (1, 65534), # 0 = NODATA, 65535 = SATURATED
            "green":    (1, 65534),
            "red":      (1, 65534),
            "nir":      (1, 65534),
            }

    masked_data = mask_invalid(data, valid_range)

    # 0 SC_NODATA # 1 SC_SATURATED_DEFECTIVE # 2 SC_DARK_FEATURE_SHADOW
    # 3 SC_CLOUD_SHADOW # 4 SC_VEGETATION # 5 SC_NOT_VEGETATED
    # 6 SC_WATER # 7 SC_UNCLASSIFIED # 8 SC_CLOUD_MEDIUM_PROBA
    # 9 SC_CLOUD_HIGH_PROBA # 10 SC_THIN_CIRRUS # 11 SC_SNOW_ICE
    pixel_qa_flags = [0, 1, 2, 3, 7, 8, 9, 10, 11]
    keep = np.invert(qa_da.isin(pixel_qa_flags))
    masked_data = masked_data.where(keep)

    scale = 1./10000. # BOA_QUANTIFICATION_VALUE
    offset = -1000 # BOA_ADD_OFFSET
    scaled_data = (masked_data + offset) * scale

    scaled_data = scaled_data.where(scaled_data <= 1.00)
    scaled_data = scaled_data.where(scaled_data >= 0.00)

    das = list()

    if np.all([x in scaled_data.band.values for x in ["red", "nir"]]):
        das.append(calc_ndvi(scaled_data))

    weights = {
        "blue": 0.074,
        "green": 0.083,
        "red": 0.334,
        "nir": 0.356,
        "offset": 0.033,
    }

    if np.all([x in scaled_data.band.values for x in ["blue", "green", "red", "nir"]]):
        das.append(calc_albedo(scaled_data, weights = weights))

    ds = xr.merge(das)

    return ds

def default_vars(product_name, req_vars):

    variables = {
        "S2MSI2A": {
                    "_B02_20m.jp2": [(), "blue"],
                    "_B03_20m.jp2": [(), "green"],
                    "_B04_20m.jp2": [(), "red"],
                    "_B8A_20m.jp2": [(), "nir"],
                    "_SCL_20m.jp2": [(), "qa"],
                },
    }

    req_dl_vars = {
        "S2MSI2A": {
            "ndvi": ["_B04_20m.jp2", "_B8A_20m.jp2", "_SCL_20m.jp2"],
            "r0": ["_B02_20m.jp2", "_B03_20m.jp2", "_B04_20m.jp2", "_B8A_20m.jp2", "_SCL_20m.jp2"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out

def default_post_processors(product_name, req_vars):
    
    post_processors = {
        "S2MSI2A": {
            "ndvi": [],
            "r0": [],
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def time_func(fn):
    dtime = np.datetime64(dt.strptime(fn.split("_")[2], "%Y%m%dT%H%M%S"))
    return dtime

def download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = None, post_processors = None, 
                extra_search_kwargs = {"cloudcoverpercentage": (0, 30)}):

    product_folder = os.path.join(folder, "SENTINEL2")

    fn = os.path.join(product_folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    if isinstance(timelim[0], str):
        timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    search_kwargs = {
                        "platformname": "Sentinel-2",
                        "producttype": product_name,
                        # "limit": 10,
    }

    search_kwargs = {**search_kwargs, **extra_search_kwargs}

    def node_filter(node_info):
        fn = os.path.split(node_info["node_path"])[-1]
        to_dl = list(variables.keys())
        return np.any([x in fn for x in to_dl])

    scenes = sentinelapi.download(product_folder, latlim, lonlim, timelim, search_kwargs, node_filter = node_filter)

    ds = sentinelapi.process_sentinel(scenes, variables, process_s2, 
                                        time_func, f"{product_name}.nc", bb = bb)

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/On My Mac/sentinel_dl_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2022-06-01", "2022-06-11"]
    product_name = 'S2MSI2A'
    req_vars = ["ndvi", "r0"]
    post_processors = None
    variables = None
    extra_search_kwargs = {"cloudcoverpercentage": (0, 30)}

    ds = download(folder, latlim, lonlim, timelim, product_name, req_vars, 
                variables = None,  post_processors = None,
                extra_search_kwargs = extra_search_kwargs
                 )
