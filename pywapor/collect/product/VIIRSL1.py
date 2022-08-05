# https://gitlab.ssec.wisc.edu/sips/OrbNavClient
# https://sips.ssec.wisc.edu/orbnav#/

import json
import fnmatch
import os
import datetime
import glob
import tqdm
import geopy.distance
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime as dt
from itertools import chain
from pywapor.general.logger import log
from pywapor.collect import accounts
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.processing_functions import save_ds, open_ds
from pywapor.general.curvilinear import create_grid, regrid
from pywapor.collect.protocol.crawler import download_url, download_urls

def regrid_VNP(workdir, latlim, lonlim, dx_dy = (0.0033, 0.0033)):

    # Reformat bounding-box.
    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    # Search for VNP02 images in workdir.
    ncs02 = glob.glob(os.path.join(workdir, "VNP02IMG.A*.nc"))

    # Make inventory of complete scenes in workdir.
    scenes = dict()
    for nc02 in ncs02:
        dt_str = ".".join(os.path.split(nc02)[-1].split(".")[1:3])
        ncs03 = glob.glob(os.path.join(workdir, f"VNP03IMG.{dt_str}*.nc"))
        if len(ncs03) != 1:
            continue
        else:
            nc03 = ncs03[0]
        dt = datetime.datetime.strptime(dt_str, "A%Y%j.%H%M")
        scenes[dt] = (nc02, nc03)

    # Create list to store xr.Datasets from a single scene.
    dss = list()

    # Loop over the scenes.
    for dt, (nc02, nc03) in tqdm.tqdm(scenes.items()):

        # Define tile output path.
        fp = os.path.join(workdir, f"VNP_{dt:%Y%j%H%M}.nc")

        if os.path.isfile(fp):
            dss.append(fp)
            continue

        # Open the datasets.
        ds1 = xr.open_dataset(nc02, group = "observation_data", decode_cf = False)
        ds2 = xr.open_dataset(nc03, group = "geolocation_data", decode_cf = False, chunks = "auto")

        # Convert DN to Brightness Temperature using the provided loopup table.
        bt_da = ds1.I05_brightness_temperature_lut.isel(number_of_LUT_values = ds1.I05)

        # Rename some things.
        ds = ds2[["latitude", "longitude"]].rename({"latitude": "y", "longitude": "x"})

        # Chunk and mask invalid pixels.
        ds["bt"] = bt_da.chunk("auto").where((bt_da >= bt_da.valid_min) & 
                                            (bt_da <= bt_da.valid_max) & 
                                            (bt_da != bt_da._FillValue))

        # Move `x` and `y` from variables to coordinates
        ds = ds.set_coords(["x", "y"])

        # Create mask for selecting the bounding-box in the data.
        buffer = 0.2
        mask = ((ds.y >= bb[1] - buffer) &
                (ds.y <= bb[3] + buffer) & 
                (ds.x >= bb[0] - buffer) & 
                (ds.x <= bb[2] + buffer))

        # Apply the mask.
        ds = ds.where(mask, drop = True)

        # Save intermediate file.
        fp_temp = os.path.join(workdir, "temp.nc")
        ds = save_ds(ds, fp_temp)

        # Create rectolinear grid.
        grid_ds = create_grid(ds, dx_dy[0], dx_dy[1], bb = bb)

        # Regrid from curvilinear to rectolinear grid.
        out = regrid(grid_ds, ds)

        # Set some metadata.
        out = out.rio.write_crs(4326)
        out = out.rio.clip_box(*bb)
        out.bt.attrs = {k: v for k, v in ds.bt.attrs.items() if k in ["long_name", "units"]}

        # Add time dimension.
        out = out.expand_dims({"time": 1}).assign_coords({"time": [dt]})

        # Save regridded tile.
        out = save_ds(out, fp)
        out = out.close()
        dss.append(fp)

        # Remove intermediate file.
        os.remove(fp_temp)

    dss = [xr.open_dataset(x, chunks = "auto", decode_coords="all") for x in dss]
    ds = xr.merge(dss)

    return ds

def boxtimes(latlim, lonlim, timelim, folder):
    """Check in which 6-min periods SUOMI NPPs swatch crosses a AOI
    defined by `latlim` and `lonlim` and remove night recordings.

    Parameters
    ----------
    latlim : list
        Latitude limits.
    lonlim : list
        Longitude limits.
    timelim : list
        Time limits.

    Returns
    -------
    list
        List containing year, doy and time at which SUOMI NPP passes over
        the bb during the day (night is filtered out).
    """

    # Convert to datetime object if necessary
    if isinstance(timelim[0], str):
        timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    # Expand bb with 1500km (is half of SUOMI NPP swath width).
    ur = (latlim[1], lonlim[1])
    ll = (latlim[0], lonlim[0])
    very_ur = geopy.distance.distance(kilometers=1500).destination(ur, bearing=45)
    very_ll = geopy.distance.distance(kilometers=1500).destination(ll, bearing=45+180)

    # Define search kwargs.
    kwargs = {
        "sat" : '37849',
        "start" : timelim[0].strftime("%Y-%m-%dT%H:%M:%SZ"),
        "end" : timelim[1].strftime("%Y-%m-%dT%H:%M:%SZ"),
        "ur" : f"{int(np.ceil(very_ur.latitude))},{int(np.ceil(very_ur.longitude))}",
        "ll" : f"{int(np.floor(very_ll.latitude))},{int(np.floor(very_ll.longitude))}",
    }

    # Search boxtimes.
    base_url = "https://sips.ssec.wisc.edu/orbnav/api/v1/boxtimes.json"
    url = base_url + "?&" + "&".join([f"{k}={v}" for k, v in kwargs.items()])
    fp = os.path.join(folder, "boxtimes.json")
    _ = download_url(url, fp)
    data = json.load(open(fp))

    # Group to SUOMI NPP 6 minute tiles.
    dates = [datetime.datetime.strptime(x[0], "%Y-%m-%dT%H:%M:%SZ") for x in chain.from_iterable(data["data"]) if x[3] > 40.]
    df = pd.DataFrame({"date": dates}).set_index("date")
    count = df.groupby(pd.Grouper(freq = "6min")).apply(lambda x: len(x))
    to_dl = count[count >= 1].index.values
    os.remove(fp)

    # Create list with relevant year/doy/times tuples.
    year_doy_time = [[pd.to_datetime(x).strftime("%Y"), pd.to_datetime(x).strftime("%j"), pd.to_datetime(x).strftime("%H%M")] for x in to_dl]

    return year_doy_time

def find_VIIRSL1_urls(year_doy_time, product, workdir, 
                        server_folders = {"VNP02IMG": 5110, "VNP03IMG": 5200}):
    """Given a list of year/doy/time tuples returns the exact urls for that
    datetime. Also check if the linked file already exists in workdir, in which
    case the url is omitted from the returned list.

    Parameters
    ----------
    year_doy_time : list
        List containing year, doy and time at which SUOMI NPP passes over
        the bb during the day (night is filtered out).
    product : str
        Name of the product to search.

    Returns
    -------
    list
        Urls linking to the product.
    """
    # Make empty url list.
    urls = list()

    # Loop over year/doy/times.
    for year, doy, time in year_doy_time:

        # Define url at which to find the json with the exact tile urls.
        url = f"https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/{server_folders[product]}/{product}/{year}/{doy}.json"

        # Download the json.
        fp = os.path.join(workdir, f"{product}_{year}{doy}.json")
        if not os.path.isfile(fp):
            _ = download_url(url, fp)

        # Open the json and find the exact url of the required scene.
        all_scenes = [x["name"] for x in json.load(open(fp))]
        req_scenes = fnmatch.filter(all_scenes, f"{product}.A{year}{doy}.{time}.*.*.nc")

        # Check if a corrct url is found.
        if len(req_scenes) == 1:
            fn = req_scenes[0]
        else:
            continue

        # Check if the file already exists in workdir, otherwise add it to `urls`.
        url = f"https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/{server_folders[product]}/{product}/{year}/{doy}/{fn}"
        fp = os.path.join(workdir, fn)
        if os.path.exists(fp):
            continue
        urls.append(url)

    return urls

def check_tiles(year_doy_time, latlim, lonlim, workdir, product = "VNP03IMG"):
    """Check if downloaded tiles actually contain data for the given AOI. If not
    removes its year/doy/time from the list.

    Parameters
    ----------
    year_doy_time : list
        List containing year, doy and time at which SUOMI NPP passes over
        the bb during the day (night is filtered out).
    latlim : list
        Latitude limits.
    lonlim : list
        Longitude limits.
    workdir : str
        Path to working directory
    product : str, optional
        name of the product to check, by default "VNP03IMG"

    Returns
    -------
    tuple
        Two lists, the first an updated `year_date_time`, the other contains the 
        values dropped from `year_date_time`.
    """
    # Create empty list to store dropped `year/doy/times`.
    dropped = list()

    # Reformat bounding-box.
    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    # Loop over `year/doy/times`.
    for year, doy, time in year_doy_time:

        # Search for the relevant .nc file.
        fps = glob.glob(os.path.join(workdir, f"{product}.A{year}{doy}.{time}.*.nc"))
        assert len(fps) == 1
        fp = fps[0]

        # Open the geolocation data.
        ds = xr.open_dataset(fp, group = "geolocation_data", chunks = "auto")

        # Mask the irrelevant pixels.
        buffer = 0.2
        mask = ((ds.latitude >= bb[1] - buffer) &
                (ds.latitude <= bb[3] + buffer) & 
                (ds.longitude >= bb[0] - buffer) & 
                (ds.longitude <= bb[2] + buffer))

        # Check if there is relevant data in this tile.
        if mask.sum().values < 100:
            year_doy_time.remove([year, doy, time])
            dropped.append([year, doy, time])

    return year_doy_time, dropped

def default_vars(product_name, req_vars):

    variables = {
        "VNP02IMG": {"bt": {(), "bt"}},
    }

    req_dl_vars = {
        "VNP02IMG": {
            "bt": ["bt"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars = None):

    post_processors = {
        "VNP02IMG": {
            "bt": []
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def download(folder, latlim, lonlim, timelim, product_name, req_vars,
                variables = None, post_processors = None):

    folder = os.path.join(folder, "VIIRSL1")
    if not os.path.exists(folder):
        os.makedirs(folder)

    fn = os.path.join(folder, f"{product_name}.nc")
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

    # Find SUOMI NPP overpass times.
    year_doy_time = boxtimes(latlim, lonlim, timelim, folder)

    # Search for valid urls at overpass times for geolocation files.
    urls = find_VIIRSL1_urls(year_doy_time, "VNP03IMG", folder)

    # Download urls.
    token, _ = accounts.get('VIIRSL1')
    headers = {'Authorization': 'Bearer ' + token}
    _ = download_urls(urls, folder, None, parallel = 3, headers = headers)

    # Check if tiles contain data in AOI.
    year_doy_time, _ = check_tiles(year_doy_time, latlim, lonlim, folder)

    # Find urls of the actual data.
    urls = find_VIIRSL1_urls(year_doy_time, product_name, folder)

    # Download the data product.
    _ = download_urls(urls, folder, None, parallel = 3, headers = headers)

    # Combine geolocations and data in rectilinear grid.
    ds = regrid_VNP(folder, latlim, lonlim, dx_dy = (0.0033, 0.0033))

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    ds = save_ds(ds, fn)

    return ds

if __name__ == "__main__":

    folder = "/Users/hmcoerver/On My Mac/viirs_test/"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2022-06-01", "2022-06-02"]
    product_name = "VNP02IMG"
    req_vars = ["bt"]

    ds = download(folder, latlim, lonlim, timelim, product_name, req_vars)