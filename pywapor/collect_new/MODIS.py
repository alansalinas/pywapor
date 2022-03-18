import datetime
import rasterio
import os
import json
import rasterio.warp
import rioxarray.merge
import xarray as xr
import numpy as np
from tqdm import tqdm
from shapely.geometry.polygon import Polygon
from pydap.cas.urs import setup_session
from shapely.geometry import shape
from pywapor.collect_new.projections import get_crss
from pywapor.collect_new.opendap import opendap_to_xarray
from pywapor.general.processing_functions import save_ds, reproject_ds

def download(folder, product_name, latlim, lonlim, timelim, un_pw):

    # Load some defaults.
    crss = get_crss()
    rename_keep_vars = default_vars(product_name)

    # Convert bounding-box to MODIS sinoidal projection.
    lonlim_source, latlim_source = rasterio.warp.transform(crss["WGS84"],
                                                            crss["MODIS"],
                                                            lonlim, latlim)
    # Determine which MODIS tiles are required.
    tiles = tiles_intersect(latlim, lonlim)

    # Create xr.Dataset selector.
    select = {"YDim": latlim_source[::-1], 
                "XDim": lonlim_source, 
                "time": [np.datetime64(timelim[0]),np.datetime64(timelim[1])]}

    dss = list()
    fps = list()

    # Loop over tiles.
    for htile, vtile in tqdm(tiles):

        # Define final output name.
        fn = f"{product_name}_h{htile:02d}v{vtile:02d}.nc"

        # Create connection with OPeNDAP tile.
        store = create_store(product_name, vtile, htile, un_pw)

        # Download the required data.
        fp_temp = os.path.join(folder, fn.replace(".nc", "_temp.nc"))
        ds_temp = opendap_to_xarray(store, fp_temp, select, rename_keep_vars)

        if isinstance(ds_temp, type(None)):
            continue

        # Reproject data to WGS84.
        fp_proj = os.path.join(folder, fn.replace(".nc", "_proj.nc"))
        dss.append(reproject_ds(ds_temp, fp_proj, crss["MODIS"], crss["WGS84"]))
        fps.append(fp_proj)

        # Remove data in original projection.
        os.remove(fp_temp)

    # Merge all the tiles.
    ds = rioxarray.merge.merge_datasets(dss)

    # Save final output.
    fp_all = os.path.join(folder, fn)
    ds = save_ds(ds, fp_all)

    # Remove the reprojected tiles.
    for fp in fps:
        os.remove(fp)

    return ds

def get_url(product, vtile, htile):
    url = f"https://opendap.cr.usgs.gov/opendap/hyrax/{product}/h{htile:02d}v{vtile:02d}.ncml"
    return url

def tiles_intersect(latlim, lonlim):
    with open(r"MODIS_tiles.geojson") as f:
        features = json.load(f)["features"]
    aoi = Polygon.from_bounds(latlim[0], lonlim[0], latlim[1], lonlim[1])
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

def create_store(product_name, vtile, htile, un_pw):
    url = get_url(product_name, vtile, htile)
    login_url = url + ".ascii?YDim[0],XDim[0],time[0]"
    session = setup_session(un_pw[0], un_pw[1], check_url=login_url)
    store = xr.backends.PydapDataStore.open(url, session=session)
    return store

def default_vars(product_name):
    rename_keep_vars = {
        "MOD13Q1.061": {"_250m_16_days_NDVI": "ndvi"},
    }
    return rename_keep_vars[product_name]

if __name__ == "__main__":

    vtile = 6
    htile = 20

    folder = r"/Users/hmcoerver/Downloads"

    product_name = "MOD13Q1.061"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 10, 11)]

    import pywapor.collect.accounts as accounts
    un_pw = accounts.get("NASA")

    ds = download(folder, product_name, latlim, lonlim, timelim, un_pw)

#%%
