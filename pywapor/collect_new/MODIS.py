"""
    t1 = datetime.datetime.now()
    dl_args = {
        "Dir": r"/Users/hmcoerver/Downloads/dl_test_old", 
        "latlim": latlim,
        "lonlim": lonlim,
        "Startdate": timelim[0],
        "Enddate": timelim[1],
        "buffer_dates": False
    }

    ndvi_files = pywapor.collect.downloader.collect_sources("ndvi", ["MOD13","MYD13"], dl_args, extra_source_locations = None)
    r0_files = pywapor.collect.downloader.collect_sources("r0", ["MCD43"], dl_args, extra_source_locations = None)
    lst_files = pywapor.collect.downloader.collect_sources("lst", ["MOD11","MYD11"], dl_args, extra_source_locations = None)
    t2 = datetime.datetime.now()
    print("old", t2 - t1)
"""
import datetime
import rasterio
import os
import json
import pywapor
import rasterio.crs
import rasterio.warp
import rioxarray.merge
import xarray as xr
import numpy as np
from functools import partial
import pywapor.collect.accounts as accounts
from shapely.geometry.polygon import Polygon
from pywapor.general.logger import log
from shapely.geometry import shape
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap
from pywapor.general.processing_functions import save_ds, create_selection
from pywapor.general import bitmasks

def download(folder, product_name, latlim, lonlim, timelim, variables = None,
                post_processors = None):

    coords = {"x": "XDim", "y": "YDim", "t": "time"}

    # Get username and password
    un_pw = accounts.get("NASA")

    # Determine which MODIS tiles are required.
    tiles = tiles_intersect(latlim, lonlim)

    # Define which variables to download.
    if isinstance(variables, type(None)):
        variables = default_vars(product_name)

    # Define which post processors to apply.
    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name)

    # Create selection object.
    selection = create_selection(latlim, lonlim, timelim, coords,
                                target_crs = get_crss("MODIS"))

    dss = list()

    for htile, vtile in tiles:

        # Define output name.
        fn = f"{product_name}_h{htile:02d}v{vtile:02d}.nc"
        fp = os.path.join(folder, fn)

        log.info(f"--> Downloading {fn}.")

        if os.path.isfile(fp):
            # Open existing dataset.
            ds_data = xr.open_dataset(fp, decode_coords = "all")
        else:
            # Define MODIS tile URL.
            base_url = get_url(product_name, htile, vtile)

            # Start OPeNDAP session.
            idxs, session = opendap.start_session(base_url, selection, un_pw)
            
            # Make data request URL.
            url_data = opendap.create_url(base_url, idxs, variables)
            
            # Download data.
            ds_data = opendap.download_url(url_data, fp, session, coords)

        # Rename (spatial) variables.
        ds_data = opendap.process_ds(ds_data, coords, variables)

        dss.append(ds_data)

    # Merge all the tiles.
    ds = rioxarray.merge.merge_datasets(dss)

    # Reproject to WGS84.
    ds = ds.rio.reproject(rasterio.crs.CRS.from_epsg(4326))

    # Apply product specific functions.
    for func in post_processors:
        ds = func(ds)

    fp = os.path.join(folder, f"{product_name}.nc")
    ds = save_ds(ds, fp, decode_coords = "all")

    return ds

def get_url(product, htile, vtile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/{product}/h{htile:02d}v{vtile:02d}.ncml.nc4?"
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
        "MOD11A1.061": [mask_bitwise_qa],#, expand_time_dim],
        "MYD11A1.061": [mask_bitwise_qa],#, expand_time_dim],
        "MCD43A3.061": [shortwave_ro, 
                        partial(mask_qa, to_mask = "r0", masker = ("r0_qa", 1.))],
    }
    return post_processors[product_name]

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/dl_test_new"

    wanted = [
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

    t1 = datetime.datetime.now()
    for product_name in wanted:
        print(product_name)
        download(folder, product_name, latlim, lonlim, timelim)
    t2 = datetime.datetime.now()
    print("new", t2 - t1)