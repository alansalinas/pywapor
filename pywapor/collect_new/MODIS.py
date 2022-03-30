import datetime
import os
import json
import pywapor
import rasterio.crs
import rasterio.warp
import numpy as np
from functools import partial
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from pywapor.general.logger import log
from shapely.geometry import shape
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap
from pywapor.general import bitmasks

def fn_func(product_name, tile):
    fn = f"{product_name}_h{tile[0]:02d}v{tile[1]:02d}.nc"
    return fn

def url_func(product_name, tile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/{product_name}/h{tile[0]:02d}v{tile[1]:02d}.ncml.nc4?"
    return url

def tiles_intersect(latlim, lonlim):
    with open(os.path.join(pywapor.collect_new.__path__[0], "MODIS_tiles.geojson")) as f:
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

def shortwave_ro(ds):
    ds["r0"] = 0.3 * ds["white_r0"] + 0.7 * ds["black_r0"]
    ds = ds.drop_vars(["white_r0", "black_r0"])
    return ds

def expand_time_dim(ds):
    groups = ds.groupby(ds.lst_hour, squeeze = True)

    def _expand_hour_dim(x):
        hour = np.timedelta64(int(x.lst_hour.median().values * 3600), "s")
        x = x.assign_coords({"hour": hour}).expand_dims("hour")
        return x

    ds_expand = groups.map(_expand_hour_dim)
    ds_expand = ds_expand.stack({"datetime": ("hour","time")})

    new_coords = [time + hour for time, hour in zip(ds_expand.time.values, ds_expand.hour.values)]

    ds_expand = ds_expand.assign_coords({"datetime": new_coords})
    ds_expand = ds_expand.rename({"datetime": "time"}).sortby("time")
    ds_expand = ds_expand.drop_vars(["lst_hour"])
    ds_expand = ds_expand.transpose("time", "y", "x")
    ds_expand = ds_expand.dropna("time", how="all")

    return ds_expand

def mask_bitwise_qa(ds, to_mask = "lst", masker = "lst_qa", 
                    product_name = "MOD11A1.061", flags = ["good_qa"]):
    flag_bits = bitmasks.MODIS_qa_translator(product_name)
    mask = bitmasks.get_mask(ds[masker].astype("uint8"), flags, flag_bits)
    ds[to_mask] = ds[to_mask].where(mask, np.nan)
    ds = ds.drop_vars([masker])
    return ds

def mask_qa(ds, to_mask = "ndvi", masker = ("ndvi_qa", 1.0)):
    ds[to_mask] = ds[to_mask].where(ds[masker[0]] != masker[1], np.nan)
    ds = ds.drop_vars(masker[0])
    return ds

def default_vars(product_name):

    vars = {

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
    return vars[product_name]

def default_post_processors(product_name):
    post_processors = {
        "MOD13Q1.061": [mask_qa],
        "MYD13Q1.061": [mask_qa],
        "MOD11A1.061": [mask_bitwise_qa, expand_time_dim],
        "MYD11A1.061": [mask_bitwise_qa, expand_time_dim],
        "MCD43A3.061": [shortwave_ro, 
                        partial(mask_qa, to_mask = "r0", masker = ("r0_qa", 1.))],
    }
    return post_processors[product_name]

def download(folder, latlim, lonlim, timelim, product_name,
                variables = None, post_processors = None):
    tiles = tiles_intersect(latlim, lonlim)
    coords = {"x": "XDim", "y": "YDim", "t": "time"}
    variables = default_vars(product_name)
    post_processors = default_post_processors(product_name)
    data_source_crs = get_crss("MODIS")
    parallel = False
    spatial_tiles = True
    un_pw = accounts.get("NASA")
    request_dims = True
    ds = opendap.download(folder, product_name, latlim, lonlim, timelim, coords, 
                variables, post_processors, fn_func, url_func, un_pw = un_pw, 
                tiles = tiles, data_source_crs = data_source_crs, parallel = parallel, 
                spatial_tiles = spatial_tiles, request_dims = request_dims)
    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/dl_test_new"

    products = [
        'MCD43A3.061',
        'MOD11A1.061',
        'MYD11A1.061',
        'MOD13Q1.061',
        'MYD13Q1.061',
    ]

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 1)]

    # MODIS.
    for product_name in products:
        ds = download(folder, latlim, lonlim, timelim, product_name)