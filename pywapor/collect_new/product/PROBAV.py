import os
from pywapor.collect import accounts
import datetime
import pandas as pd
import pywapor.collect.accounts as accounts
from pywapor.general.logger import log
import re
from functools import partial
import requests
from pywapor.collect_new.protocol.requests import download_urls, crawl
import pywapor.general.bitmasks as bm
import xarray as xr
import numpy as np
from pywapor.general.processing_functions import process_ds, save_ds, open_ds
import tqdm
import datetime
import rioxarray.merge

def download(folder, latlim, lonlim, timelim, product_name,
                variables = None, post_processors = None):

    dates = pd.date_range(timelim[0], timelim[1], freq="D")

    bb = (lonlim[0], latlim[0], lonlim[1], latlim[1])

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name)

    if isinstance(variables, type(None)):
        variables = default_vars(product_name)

    coords = {"x": ["lon", None], "y": ["lat", None]}

    base_url = f"https://www.vito-eodata.be/PDF/datapool/Free_Data/PROBA-V_100m/{product_name}"
    coord_req = "?coord=" + ",".join([str(x) for x in bb])
    url = os.path.join(base_url, coord_req)

    session = requests.sessions.Session()
    session.auth = accounts.get("VITO")

    # Scrape urls.
    urls = find_tiles(url, dates, session)

    # Download .HDF5 tiles.
    fps = download_urls(urls, folder, session, parallel = 4)

    # Convert to netcdf.
    dss = dict()
    for fp in tqdm.tqdm(fps):
        date = datetime.datetime.strptime(fp.split("_")[-3], "%Y%m%d")
        ds = open_hdf5_groups(fp, variables, coords).expand_dims({"time": [date]})
        if date in dss.keys():
            dss[date].append(ds)
        else:
            dss[date] = [ds]

    # Merge tiles.
    dss0 = list()
    for date, datasets in dss.items():
        log.info(f"--> Merging {len(datasets)} tiles for {date.strftime('%Y-%m-%d')}.")
        bb = (lonlim[0], latlim[0], lonlim[1], latlim[1])
        ds = rioxarray.merge.merge_datasets(datasets)
        ds = ds.rio.clip_box(*bb) # NOTE using this seperately, because `bounds`-kw for `merge_datasets` bugs.
        dss0.append(ds)

    # Merge dates.
    ds = xr.concat(dss0, dim = "time")

    # Apply product specific functions.
    for func in post_processors:
        ds = func(ds)

    ds = save_ds(ds, os.path.join(folder, f"{product_name}.nc"), decode_coords = "all")

    for fp in fps:
        os.remove(fp.replace(".HDF5", ".nc"))

    return ds

def open_hdf5_groups(fp, variables, coords):

    nc_fp = fp.replace(".HDF5", ".nc")
    if os.path.isfile(nc_fp):
        ds = open_ds(nc_fp, "all")
    else:
        log.info(f"--> Converting {os.path.split(fp)[-1]} to netcdf.")
        ds = xr.open_dataset(fp, chunks = "auto")

        spatial_ref_name = [k for k, v in variables.items() if v[1] == "spatial_ref"][0]
        if (spatial_ref_name in list(ds.variables)) and (spatial_ref_name not in list(ds.coords)):
            ds = ds.set_coords((spatial_ref_name))

        for k in variables.keys():
            if k in ds.variables:
                continue
            with xr.open_dataset(fp, group = k, engine = "netcdf4", chunks = "auto") as ds_grp:
                vrs = list(ds_grp.variables)[0] # NOTE assuming each group has only 1 variable!
                ds[k] = ds_grp.rename({vrs: k})[k]

        ds = process_ds(ds, coords, variables)

        ds = save_ds(ds, nc_fp, decode_coords = "all")

    return ds

def mask_bitwise_qa(ds, to_mask, flags):
    flag_bits = bm.PROBAV_qa_translator()
    mask = bm.get_mask(ds["qa"].astype("uint8"), flags, flag_bits)
    ds[to_mask] = ds[to_mask].where(~mask, np.nan)
    return ds

def calc_r0(ds):
    ds["r0"] = 0.429 * ds["blue"] + 0.333 * ds["red"] + 0.133 * ds["nir"] + 0.105 * ds["swir"]
    return ds

def drop_vars(ds, to_drop):
    ds = ds.drop_vars(to_drop)
    return ds

def default_post_processors(product_name):
    post_processors = {
        "S5_TOC_100_m_C1": [calc_r0, 
                            partial(mask_bitwise_qa, to_mask = "r0", 
                            flags = ["bad BLUE", "bad RED", "bad NIR", "bad SWIR", 
                                    "sea", "undefined", "cloud", "ice/snow", "shadow"]),
                            partial(mask_bitwise_qa, to_mask = "ndvi", 
                            flags = ["bad RED", "bad NIR", "sea", "undefined", 
                            "cloud", "ice/snow", "shadow"]),
                            partial(drop_vars, to_drop = ["blue", "nir", "red", 
                            "swir", "vnir_vza", "swir_vza", "qa"])
                            ]
    }
    return post_processors[product_name]

def default_vars(product_name):
    vars = {
        "S5_TOC_100_m_C1": {
            "LEVEL3/NDVI":              [("lon", "lat"), "ndvi"],
            "LEVEL3/RADIOMETRY/BLUE":   [("lon", "lat"), "blue"],
            "LEVEL3/RADIOMETRY/NIR":    [("lon", "lat"), "nir"],
            "LEVEL3/RADIOMETRY/RED":    [("lon", "lat"), "red"],
            "LEVEL3/RADIOMETRY/SWIR":   [("lon", "lat"), "swir"],
            "LEVEL3/GEOMETRY/VNIR":     [("lon", "lat"), "vnir_vza"],
            "LEVEL3/GEOMETRY/SWIR":     [("lon", "lat"), "swir_vza"],
            "LEVEL3/QUALITY":           [("lon", "lat"), "qa"],
            "crs":                      [(), "spatial_ref"]
        },
    }
    return vars[product_name]

def find_tiles(url, dates, session):

    log.info("--> Searching PROBAV tiles.")

    regex = "https:.*\/\d{4}\/\?coord="
    filter_regex = "\d{4}(?=\/\?coord=)"
    urls = {"_": url}
    label_filter = dates.strftime("%Y")
    years = crawl(urls, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\?coord="
    filter_regex = "\d{4}\/\d{2}(?=\/\?coord=)"
    label_filter = dates.strftime("%Y%m")
    months = crawl(years, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\d{2}\/\?coord="
    filter_regex = "\d{4}\/\d{2}\/\d{2}(?=\/\?coord=)"
    label_filter = dates.strftime("%Y%m%d")
    days = crawl(months, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\d{2}\/.*\/\?coord="
    filter_regex = "\d{4}\/\d{2}\/\d{2}(?=\/.*\/\?coord=)"
    label_filter = dates.strftime("%Y%m%d")
    prods = crawl(days, regex, filter_regex, session, label_filter = label_filter)

    regex = ".*\.HDF5"
    filter_regex = "\d{8}(?=.*\.HDF5)"
    fns = crawl(prods, regex, filter_regex, session, list_out = True)

    dl_urls = [os.path.join(re.sub("\?coord=.*","",prods[date_str]), fn) for date_str, fn, in fns]

    log.info(f"--> Found {len(dl_urls)} PROBAV tiles.")

    return dl_urls

if __name__ == "__main__":

    product_name = "S5_TOC_100_m_C1"

    folder = r"/Users/hmcoerver/Downloads/merra2"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    timelim = [datetime.date(2020, 7, 20), datetime.date(2020, 7, 22)]

    variables = None
    post_processors = None

    ds = download(folder, latlim, lonlim, timelim, product_name,
                variables = variables, post_processors = post_processors)
