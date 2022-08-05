import datetime
import os
import json
import pywapor
import numpy as np
from functools import partial
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from pywapor.general.logger import log
from shapely.geometry import shape
from pywapor.collect.protocol.projections import get_crss
import pywapor.collect.protocol.opendap as opendap
from pywapor.general.processing_functions import open_ds
from pywapor.general import bitmasks
import pandas as pd
import warnings

def fn_func(product_name, tile):
    fn = f"{product_name}_h{tile[0]:02d}v{tile[1]:02d}.nc"
    return fn

def url_func(product_name, tile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/{product_name}/h{tile[0]:02d}v{tile[1]:02d}.ncml.nc4?"
    return url

def tiles_intersect(latlim, lonlim):
    with open(os.path.join(pywapor.collect.__path__[0], "product/MODIS_tiles.geojson")) as f:
        features = json.load(f)["features"]
    aoi = Polygon.from_bounds(lonlim[0], latlim[0], lonlim[1], latlim[1])
    tiles = list()
    for feature in features:
        shp = shape(feature["geometry"])
        tile = feature["properties"]["Name"]
        if shp.intersects(aoi):
            h, v = tile.split(" ")
            htile = int(h.split(":")[-1])
            vtile = int(v.split(":")[-1])
            tiles.append((htile, vtile))
    return tiles

def shortwave_r0(ds, *args):
    ds["r0"] = 0.3 * ds["white_r0"] + 0.7 * ds["black_r0"]
    ds = ds.drop_vars(["white_r0", "black_r0"])
    return ds

def expand_time_dim(ds, *args):

    groups = ds.groupby(ds.lst_hour, squeeze = True)

    def _expand_hour_dim(x):
        hour = np.timedelta64(int(np.nanmedian(x.lst_hour.values) * 3600), "s")
        x = x.assign_coords({"hour": hour}).expand_dims("hour")
        return x

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Slicing with an out-of-order index")

        ds_expand = groups.map(_expand_hour_dim)

        ds_expand = ds_expand.stack({"datetime": ("hour","time")})

        new_coords = [time + hour for time, hour in zip(ds_expand.time.values, ds_expand.hour.values)]
        
        try: # new versions of xarray require to drop all dimensions of a multi-index
            ds_expand = ds_expand.drop_vars(["datetime", "hour", "time"])
        except ValueError: # old versions throw an error when trying to drop sub-dimensions of a multiindex.
            ds_expand = ds_expand.drop_vars(["datetime"])
        
        ds_expand = ds_expand.assign_coords({"datetime": new_coords}).rename({"datetime": "time"}).sortby("time")
        ds_expand = ds_expand.drop_vars(["lst_hour"])
        ds_expand = ds_expand.transpose("time", "y", "x")
        ds_expand = ds_expand.dropna("time", how="all")

    return ds_expand

def mask_bitwise_qa(ds, var, masker = "lst_qa", 
                product_name = "MOD11A1.061", flags = ["good_qa"]):

    new_data = ds[var]

    flag_bits = bitmasks.MODIS_qa_translator(product_name)
    mask = bitmasks.get_mask(ds[masker].astype("uint8"), flags, flag_bits)
    new_data = ds[var].where(mask, np.nan)
    ds = ds.drop_vars([masker])

    ds[var] = new_data

    return ds

def mask_qa(ds, var, masker = ("ndvi_qa", 1.0)):

    new_data = ds[var].where(ds[masker[0]] != masker[1], np.nan)
    ds = ds.drop_vars(masker[0])

    ds[var] = new_data

    return ds

def default_vars(product_name, req_vars):

    variables = {

        "MOD13Q1.061": {
                    "_250m_16_days_NDVI": [("time", "YDim", "XDim"), "ndvi"],
                    "_250m_16_days_pixel_reliability": [("time", "YDim", "XDim"), "ndvi_qa"],
                    "MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection": [(), "spatial_ref"],
                        },
        "MYD13Q1.061": {
                    "_250m_16_days_NDVI": [("time", "YDim", "XDim"), "ndvi"],
                    "_250m_16_days_pixel_reliability": [("time", "YDim", "XDim"), "ndvi_qa"],
                    "MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection": [(), "spatial_ref"],
                        },
        "MOD11A1.061": {
                    "LST_Day_1km": [("time", "YDim", "XDim"), "lst"],
                    "Day_view_time": [("time", "YDim", "XDim"), "lst_hour"],
                    "QC_Day": [("time", "YDim", "XDim"), "lst_qa"],
                    "MODIS_Grid_Daily_1km_LST_eos_cf_projection": [(), "spatial_ref"],
                        },
        "MYD11A1.061": {
                    "LST_Day_1km": [("time", "YDim", "XDim"), "lst"],
                    "Day_view_time": [("time", "YDim", "XDim"), "lst_hour"],
                    "QC_Day": [("time", "YDim", "XDim"), "lst_qa"],
                    "MODIS_Grid_Daily_1km_LST_eos_cf_projection": [(), "spatial_ref"],
                        },
        "MCD43A3.061": {
                    "Albedo_WSA_shortwave": [("time", "YDim", "XDim"), "white_r0"],
                    "Albedo_BSA_shortwave": [("time", "YDim", "XDim"), "black_r0"],
                    "BRDF_Albedo_Band_Mandatory_Quality_shortwave": [("time", "YDim", "XDim"), "r0_qa"],
                    "MOD_Grid_BRDF_eos_cf_projection": [(), "spatial_ref"]
        }
    }

    req_dl_vars = {
        "MOD13Q1.061": {
            "ndvi": ["_250m_16_days_NDVI", "_250m_16_days_pixel_reliability", "MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection"],
        },
        "MYD13Q1.061": {
            "ndvi": ["_250m_16_days_NDVI", "_250m_16_days_pixel_reliability", "MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection"],
        },
        "MOD11A1.061": {
            "lst": ["LST_Day_1km", "Day_view_time", "QC_Day", "MODIS_Grid_Daily_1km_LST_eos_cf_projection"],
        },
        "MYD11A1.061": {
            "lst": ["LST_Day_1km", "Day_view_time", "QC_Day", "MODIS_Grid_Daily_1km_LST_eos_cf_projection"],
        },
        "MCD43A3.061": {
            "r0": ["Albedo_WSA_shortwave", "Albedo_BSA_shortwave", "BRDF_Albedo_Band_Mandatory_Quality_shortwave", "MOD_Grid_BRDF_eos_cf_projection"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars = None):

    post_processors = {
        "MOD13Q1.061": {
            "ndvi": [mask_qa]
            },
        "MYD13Q1.061": {
            "ndvi": [mask_qa]
            },
        "MOD11A1.061": {
            "lst": [mask_bitwise_qa, expand_time_dim]
            },
        "MYD11A1.061": {
            "lst": [mask_bitwise_qa, expand_time_dim]
            },
        "MCD43A3.061": {
            "r0": [
                    shortwave_r0, 
                    partial(mask_qa, masker = ("r0_qa", 1.)),
                    ]
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def download(folder, latlim, lonlim, timelim, product_name, req_vars,
                variables = None, post_processors = None):
    folder = os.path.join(folder, "MODIS")

    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    if product_name == "MOD13Q1.061" or product_name == "MYD13Q1.061":
        timedelta = np.timedelta64(8, "D")
        timelim[0] = timelim[0] - pd.Timedelta(timedelta)
    elif product_name == "MCD43A3.061":
        timedelta = np.timedelta64(12, "h")
        timelim[0] = timelim[0] - pd.Timedelta(timedelta)
    else:
        timedelta = None

    tiles = tiles_intersect(latlim, lonlim)
    coords = {"x": ["XDim", lonlim], "y": ["YDim", latlim], "t": ["time",timelim]}
    
    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    data_source_crs = get_crss("MODIS")
    parallel = False
    spatial_tiles = True
    un_pw = accounts.get("NASA")
    request_dims = True
    ds = opendap.download(folder, product_name, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims,
                timedelta = timedelta)
    return ds

if __name__ == "__main__":

    products = [
        # ('MCD43A3.061', ["r0"]),
        ('MOD11A1.061', ["lst"]),
        # ('MYD11A1.061', ["lst"]),
        # ('MOD13Q1.061', ["ndvi"]),
        # ('MYD13Q1.061', ["ndvi"]),
    ]

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 11)]


    for product_name, req_vars in products:
        ds = download(folder, latlim, lonlim, timelim, product_name, req_vars)
        print(ds.rio.crs, ds.rio.grid_mapping)