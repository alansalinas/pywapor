import os
import pywapor
import cdsapi
import logging
import pandas as pd
import numpy as np
import shutil
import glob
import xarray as xr
import rasterio
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

    # Loop over requests
    for setting in settings:
        ext = {"zip": "zip", "netcdf": "nc", "grib": "grib", "tgz": "tar.gz"}[setting["format"]]
        fn = f"{setting['year']}_{setting['month']}_{setting['variable']}"
        subfolder = os.path.join(folder, fn)
        fp = os.path.join(subfolder, f"{fn}.{ext}")

        if not os.path.exists(subfolder):
            os.makedirs(subfolder)

        # Make the request
        _ = c.retrieve(product_name, setting, fp)

        # Unpack if necessary
        if ext in ["zip", "tar.gz"]:
            shutil.unpack_archive(fp, subfolder)
            os.remove(fp)

        # Rename to common variable name in xr.Dataset
        ncs = glob.glob(os.path.join(subfolder, "*.nc"))
        ds = xr.open_mfdataset(ncs)
        ds = ds.rename_vars({list(ds.data_vars)[0]: variables[setting["variable"]][1]})
        dss.append(ds)
            
    # Merge everything together.
    ds = xr.merge(dss)

    # Clean up the dataset.
    ds = ds.rename_dims({"lat": "y", "lon": "x"}).rename_vars({"lat": "y", "lon": "x"})
    ds = ds.rio.write_crs(rasterio.crs.CRS.from_epsg(4326))
    ds = ds.rio.write_grid_mapping("spatial_ref")
    for var in list(ds.data_vars):
        ds[var].attrs = {k:v for k,v in ds[var].attrs.items() if k == "units"}
    ds = ds.sortby("y", ascending = False)
    ds = ds.sortby("x")
    ds.attrs = {}

    # Save the netcdf.
    ds = save_ds(ds, fn_final)
    
    return ds


