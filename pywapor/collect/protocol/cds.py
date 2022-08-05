import os
import pywapor
import cdsapi
import logging
import pandas as pd
import numpy as np
import glob
import xarray as xr
import rasterio.crs
import re
import shutil
from pywapor.general.processing_functions import save_ds

def create_time_settings(timelim):
    dates = pd.date_range(timelim[0], timelim[1], freq = "D")
    settings = list()
    for yr in np.unique(dates.year):
        for mnth in np.unique(dates.month[dates.year == yr]):
            days = dates.day[np.all([dates.year == yr, dates.month == mnth], axis = 0)]
            settings.append({"year": f"{yr}", "month": f"{mnth:02d}", "day": [f"{x:02d}" for x in days]})
    return settings

def download(folder, product_name, latlim, lonlim, timelim, variables):

    fn_final = os.path.join(folder, f"{product_name}.nc")

    # Create the settings for each individual request.
    time_settings = create_time_settings(timelim)
    area_settings = {"area": [latlim[1], lonlim[0], latlim[0], lonlim[1]]}
    settings = list()
    for var, extra_settings in variables.items():
        for t_setting in time_settings:
            settings.append({**t_setting, **extra_settings[0], 
                                **{"variable": var}, **area_settings})

    # Load api key.
    url, key = pywapor.collect.accounts.get("ECMWF")

    # Connect to server.
    c = cdsapi.Client(url = url, key = key, verify = True)
    log_settings = logging.getLogger("cdsapi")
    log_settings.setLevel(logging.CRITICAL)

    dss = list()
    subfolders = list()

    # Loop over requests
    for setting in settings:

        ext = {"zip": "zip", "netcdf": "nc", "grib": "grib", "tgz": "tar.gz"}[setting["format"]]
        fn = f"{setting['year']}_{setting['month']}_{setting['variable']}_{product_name}"
        fp = os.path.join(folder, f"{fn}.{ext}")

        # Make the request
        if not os.path.isfile(fp):
            _ = c.retrieve(product_name, setting, fp)

        # Unpack if necessary
        if ext in ["zip", "tar.gz"]:
            subfolder = os.path.join(folder, fn)
            if not os.path.exists(subfolder):
                os.makedirs(subfolder)
            shutil.unpack_archive(fp, subfolder)
            fps = glob.glob(os.path.join(subfolder, "*.nc"))
            subfolders.append(subfolder)
        else:
            fps = [fp]

        # Open downloaded data
        ds = xr.open_mfdataset(fps)

        das = list()

        time_offset = {"sis-agrometeorological-indicators": 12,
                "reanalysis-era5-single-levels": 0}[product_name]

        for var in ds.data_vars:
            # Fix time of relative humidity in agERA5.
            if bool(re.search(r'_[01]\dh', var)):
                offset = int(re.search(r'_[01]\dh', var).group()[1:-1])
                da = ds[var].assign_coords({"time": ds[var].time + np.timedelta64(offset, "h")})
            # Adjust time to middle of day for daily data.
            else:
                da = ds[var].assign_coords({"time": ds[var].time + np.timedelta64(time_offset, "h")})
            das.append(da)

        ds = xr.concat(das, dim="time").to_dataset().sortby("time")

        renames = {x: variables[setting["variable"]][1] for x in ds.data_vars}
        ds = ds.rename_vars(renames)
        
        dss.append(ds)

    # Merge everything together.
    ds = xr.merge(dss)

    # Clean up the dataset.
    relevant_coords = {
        "lat": "y", 
        "latitude": "y", 
        "lon": "x", 
        "longitude": "x", 
        # "time": "time",
    }

    coord_renames = {k: v for k, v in relevant_coords.items() if k in ds.coords}
    ds = ds.rename_dims(coord_renames)
    ds = ds.rename_vars(coord_renames)
    ds = ds.drop_vars([x for x in ds.coords if x not in ds.dims])
    ds = ds.rio.write_crs(rasterio.crs.CRS.from_epsg(4326))
    ds = ds.rio.write_grid_mapping("spatial_ref")
    for var in list(ds.data_vars):
        ds[var].attrs = {k:v for k,v in ds[var].attrs.items() if k == "units"}
    ds = ds.sortby("y", ascending = False)
    ds = ds.sortby("x")
    ds.attrs = {}

    # Save the netcdf.
    ds = save_ds(ds, fn_final)

    # Remove unpacked zips.
    for subfolder in subfolders:
        shutil.rmtree(subfolder)
    
    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/On My Mac/era_test/ERA5"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2021-06-26", "2021-07-11"]

    # product_name = "sis-agrometeorological-indicators"
    product_name = "reanalysis-era5-single-levels"

    # req_vars = ["t_air", "t_dew", "rh", "u", "vp", "ra"]
    req_vars = ["u_10m", "v_10m", "t_dew", "p_air_0", "p_air", "t_air"]

    variables = pywapor.collect.product.ERA5.default_vars(product_name, req_vars)
