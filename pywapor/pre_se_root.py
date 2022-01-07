import tqdm
import warnings
import os
import copy
import xarray as xr
import numpy as np
import pandas as pd
import pywapor.collect as c
import pywapor.general as g
from pywapor.collect.downloader import collect_sources
from pywapor.general.logger import log, adjust_logger
from datetime import datetime as dat
from datetime import time as datt
from pywapor.general.compositer import check_geots, preprocess_func
from pywapor.enhancers.temperature import kelvin_to_celsius, lapse_rate
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.compositer import calculate_ds
from pywapor.general import processing_functions as pf

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1",
        extra_sources = None, extra_source_locations = None):

    log_write = True
    log_level = "INFO"
    adjust_logger(log_write, project_folder, log_level)

    log.info("> PRE_SE_ROOT").add()

    # Disable Warnings
    warnings.filterwarnings('ignore')

    # Load required variable sources.
    levels = g.variables.get_source_level_selections()

    if isinstance(level, str):
        source_selection = levels[level]
        level_name = level
    elif isinstance(level, dict):
        source_selection = copy.copy(level)
        if "level_name" in source_selection.keys():
            level_name = source_selection.pop("level_name")
        else:
            level_name = "custom"

    raw_folder = os.path.join(project_folder, "RAW")
    level_folder = os.path.join(project_folder, level_name)
    temp_folder = os.path.join(project_folder, "temporary")

    dl_args = {
                "Dir": raw_folder, 
                "latlim": latlim, 
                "lonlim": lonlim, 
                "Startdate": startdate, 
                "Enddate": enddate,
                }

    #### NDVI ####
    log.info(f"# ndvi")
    raw_ndvi_files = collect_sources("ndvi", source_selection["ndvi"], dl_args, extra_source_locations)

    cmeta = {
        "composite_type": False,
        "temporal_interp": False,
        "temporal_interp_freq": 1,
        "spatial_interp": "nearest",
        "var_name": "ndvi",
        "var_unit": "-",
    }
    ds_ndvi = g.compositer.main(cmeta, raw_ndvi_files, None, temp_folder, None, lean_output = False)

    # example_ds = ds_ndvi.isel(time = 0).drop_vars(["time", "sources"])
    example_fh, example_ds, example_geoinfo, resolution = pf.select_template(raw_ndvi_files)

    #### LST ####
    log.info("# lst")
    raw_lst_files = collect_sources("lst", source_selection["lst"], dl_args, extra_source_locations)

    ds_lst = combine_lst(raw_lst_files)

    #### METEO ####
    log.info("> METEO").add()

    all_vars = ["t_air_i", "u2m_i", "v2m_i", "qv_i", 
                "wv_i", "p_air_i", "p_air_0_i"]

    dss = list()

    meteo_enhancements = {
        ("MERRA2",  "t_air_i"):        [kelvin_to_celsius],
        ("GEOS5",   "t_air_i"):        [kelvin_to_celsius],
    }

    for var in all_vars:

        meteo_source = source_selection[var][0]
        meteo_freq = {"GEOS5": "3H", "MERRA2": "H"}[meteo_source]

        sds, eds, prds = calc_periods(ds_lst.time.values, meteo_freq)

        dl_args["Startdate"] = sds
        dl_args["Enddate"] = eds
        dl_args["Periods"] = prds

        db = collect_sources(var, [meteo_source], dl_args, 
                                extra_source_locations = extra_source_locations)[0]

        check_geots(db)
        ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", 
                                combine = "nested", preprocess = preprocess_func)
        ds = ds.rename({"band_data": var})

        if (meteo_source, var) in meteo_enhancements.keys():
            for enhancer in meteo_enhancements[(meteo_source, var)]:
                ds, label = apply_enhancer(ds, var, enhancer, 
                                                source = meteo_source)
                ds, _ = calculate_ds(ds, label = label)

        dss.append(ds)

    ds_meteo = xr.merge(dss, combine_attrs = "drop")

    log.sub().info("< METEO")

    # spatial interpolation
    ds_lst2 = ds_lst.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
    ds_meteo2 = ds_meteo.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
    ds_ndvi2 = ds_ndvi.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)

    # temporal interpolation
    ds_meteo3 = ds_meteo2.interp(time = ds_lst2.time, method = "nearest", kwargs={"fill_value": "extrapolate"},) # TODO download more meteo data and switch to linear
    ds_ndvi3 = ds_ndvi2.interp(time = ds_lst.time, method = "linear", kwargs={"fill_value": "extrapolate"},)

    ds_se_root = xr.merge([ds_lst2, ds_meteo3, ds_ndvi3])

    ds_se_root = g.variables.fill_attrs(ds_se_root)

    ds_se_root.attrs["geotransform"] = example_geoinfo[0]
    ds_se_root.attrs["projection"] = example_geoinfo[1]
    ds_se_root.attrs["pixel_size"] = resolution
    ds_se_root.attrs["example_file"] = example_fh

    fh = os.path.join(level_folder, "se_root_input.nc")
    ds, fh = calculate_ds(ds_se_root, fh, "--> Resampling datasets.")

    log.sub().info("< PRE_SE_ROOT")

    return ds, fh

def calc_periods(times, freq):

    sds = list()
    eds = list()
    prds = list()

    for now_time in times:

        date = pd.Timestamp(now_time).date()
        starttime = dat(date.year, date.month, date.day, 0, 0)
        endtime = dat(date.year, date.month, date.day, 23, 59)

        offsets = {"3H": 90, "H": 30}

        DateTime = pd.date_range(starttime, endtime, freq=freq) + pd.offsets.Minute(offsets[freq])
        Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - now_time))
        period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1

        sds.append(starttime)
        eds.append(endtime)
        prds.append([period])

    return sds, eds, prds

def combine_lst(lst_files):
    ds = None

    for lst_db in lst_files:
        for date, files in lst_db.items():

            ds_time = xr.open_dataset(files[1]).squeeze("band").rename({"band_data": "time", "x": "lon", "y": "lat"})
            ds_angle = xr.open_dataset(files[2]).squeeze("band").rename({"band_data": "angle", "x": "lon", "y": "lat"})
            ds_lst = xr.open_dataset(files[0]).squeeze("band").rename({"band_data": "lst", "x": "lon", "y": "lat"})

            unique_times = np.unique(ds_time.time)
            unique_times = unique_times[np.isfinite(unique_times)]

            offset_GTM = int(round(np.median(ds_time.lon) * 24 / 360)) # TODO Make this spatial for areas spanning multiple timezones

            for time in unique_times:

                hours = np.floor(time).astype(int) - offset_GTM
                minutes = np.round((time - offset_GTM - hours) * 60).astype(int)
                dt = dat.combine(date, datt(hours, minutes))

                da_angle = ds_angle.angle.where(ds_time.time == time)
                da_angle = da_angle.assign_coords({"time": dt}).expand_dims("time", axis = 0)

                da_lst = ds_lst.lst.where(ds_time.time == time, np.nan)
                da_lst = da_lst.assign_coords({"time": dt}).expand_dims("time", axis = 0)

                if isinstance(ds, type(None)):
                    ds = xr.Dataset({"angle": da_angle, "lst": da_lst})
                else:
                    ds_temp = xr.Dataset({"angle": da_angle, "lst": da_lst})
                    ds = xr.concat([ds, ds_temp], dim = "time")

    ds.attrs["time"] = "<UTC> time"

    return ds

# if __name__ == "__main__":

    # project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    # startdate = "2021-07-01"
    # enddate = "2021-07-11"
    # composite_length = "DEKAD"
    # level = "level_1"
    # extra_sources = None
    # extra_source_locations = None

    # diagnostics = { # label          # lat      # lon
    #                 "water":	    (29.44977,	30.58215),
    #                 "desert":	    (29.12343,	30.51222),
    #                 "agriculture":	(29.32301,	30.77599),
    #                 "urban":	    (29.30962,	30.84109),
    #                 }
    # # diagnostics = None

    # level = {'METEO': ['MERRA2'],
    #         'NDVI': ['MYD13'],
    #         'ALBEDO': ['MCD43'],
    #         'LST': ['MOD11', 'MYD11'],
    #         'LULC': ['WAPOR'],
    #         'DEM': ['SRTM'],
    #         'PRECIPITATION': ['CHIRPS'],
    #         'SOLAR_RADIATION': ['MERRA2'],
    #         "name": "level_1"}

    # main(project_folder, startdate, enddate, latlim, lonlim, level = level,
    #     extra_sources = extra_sources, extra_source_locations = extra_source_locations)