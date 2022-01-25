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
from pywapor.enhancers.temperature import kelvin_to_celsius
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.compositer import calculate_ds
from pywapor.general import processing_functions as pf

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1", 
            extra_source_locations = None):

#%%

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

    example_fh, example_ds, example_geoinfo, resolution = pf.select_template(raw_ndvi_files)

    cmeta = {
        "composite_type": None, # 
        "temporal_interp": False,
        "temporal_interp_freq": 1,
        "spatial_interp": "nearest",
        "var_name": "ndvi",
        "var_unit": "-",
    }
    ds_ndvi = g.compositer.main(cmeta, raw_ndvi_files, None, temp_folder, None, lean_output = False)
    ds_ndvi = ds_ndvi.drop_vars(["sources"]).chunk("auto")

#%%
    #### LST ####
    log.info("# lst")
    raw_lst_files = collect_sources("lst", source_selection["lst"], dl_args, extra_source_locations)

    # TODO `combine_lst` needs to go into collect, then this whole block can go away!
    raw_modis_files = [x for x in raw_lst_files if isinstance(x, dict)]
    raw_nonmodis_files = [x for x in raw_lst_files if isinstance(x, list)]

    if len(raw_modis_files) > 0:
        ds_lst = combine_lst(raw_lst_files) 

        ds_lst2a = ds_lst.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
        fh = os.path.join(temp_folder, "lst_ts.nc")
        ds_lst2a, _ = calculate_ds(ds_lst2a, fh, label = "--> Resampling datasets.")
    else:
        ds_lst2a = False

    if len(raw_nonmodis_files) > 0:
        cmeta = {
            "composite_type": None, # 
            "temporal_interp": False,
            "temporal_interp_freq": 1,
            "spatial_interp": "nearest",
            "var_name": "lst",
            "var_unit": "-",
        }
        ds_lst2b = g.compositer.main(cmeta, raw_nonmodis_files, None, temp_folder, example_ds, lean_output = False)
        ds_lst2b = ds_lst2b.drop_vars(["sources"]).chunk("auto")
    else:
        ds_lst2b = False

    if ds_lst2a and ds_lst2b:
        ds_lst2 = xr.merge([ds_lst2a, ds_lst2b])
    elif ds_lst2a:
        ds_lst2 = ds_lst2a
    elif ds_lst2b:
        ds_lst2 = ds_lst2b
    else:
        raise ValueError

#%%
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

        log.info(f"# {var}")

        meteo_source = source_selection[var][0]

        if len(source_selection[var]) > 1:
            log.warning(f"! --> Multiple sources not supported for '{var}', using '{meteo_source}' only.")

        meteo_freq = {"GEOS5": "3H", "MERRA2": "H"}[meteo_source]

        sds, eds, prds = calc_periods(ds_lst2.time.values, meteo_freq)

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
        else:
            ds, _ = calculate_ds(ds, fh = os.path.join(temp_folder, f"{var}_ts.nc"))

        ds = ds.chunk("auto")

        attributes = {
                "sources": [meteo_source],
            }
        ds[var].attrs = attributes

        # temporal
        ds, x, y = g.compositer.temporal_interpolation(ds, ds_lst2.time.values, "nearest")
        ds = ds.sel(time = ds_lst2.time.values)

        # spatial
        ds = ds.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
        
        dss.append(ds)

    ds_meteo3 = xr.merge(dss, combine_attrs = "drop_conflicts")

    log.sub().info("< METEO")

    ds_ndvi3, x, y = g.compositer.temporal_interpolation(ds_ndvi, ds_lst2.time.values, "linear")
    ds_ndvi3 = ds_ndvi3.sel(time = ds_lst2.time)

    ds_se_root = xr.merge([ds_lst2, ds_meteo3, ds_ndvi3])

    if "spatial_ref" in list(ds_se_root.coords):
        ds_se_root = ds_se_root.drop_vars("spatial_ref")
    if "band" in list(ds_se_root.coords):
        ds_se_root = ds_se_root.drop_vars("band")

    ds_se_root = g.variables.fill_attrs(ds_se_root)

    ds_se_root.attrs["geotransform"] = example_geoinfo[0]
    ds_se_root.attrs["projection"] = example_geoinfo[1]
    ds_se_root.attrs["pixel_size"] = resolution
    ds_se_root.attrs["example_file"] = example_fh

    fh = os.path.join(level_folder, "se_root_input.nc")
    ds, fh = calculate_ds(ds_se_root, fh, "--> Interpolating datasets.")

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

    sources = list()

    for lst_db in lst_files:

        fn = os.path.split(list(lst_db.values())[0][0])[-1]
        source = fn.split("_")[1]
        sources.append(source)

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

    ds["lst"].attrs = {"sources": sources}
    ds["angle"].attrs = {"sources": sources}

    return ds

if __name__ == "__main__":

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"
    composite_length = "DEKAD"
    # level = "level_2"
    # extra_source_locations = None

    level = {
        # Main inputs
        "ndvi":         ["LS7"],
        "r0":           ["LS7"],
        "lst":          ["LS7"],
        "lulc":         ["WAPOR"],
        "z":            ["SRTM"],
        "p_24":         ["CHIRPS"],
        "ra_24":        ["MERRA2"],

        # Daily meteo 
        't_air_24':     ["GEOS5"],
        't_air_min_24': ["GEOS5"], 
        't_air_max_24': ["GEOS5"],
        'u2m_24':       ["GEOS5"],
        'v2m_24':       ["GEOS5"],
        'p_air_0_24':   ["GEOS5"],
        'qv_24':        ["GEOS5"],

        # Instanteneous meteo
        "t_air_i":      ["GEOS5"],
        "u2m_i":        ["GEOS5"],
        "v2m_i":        ["GEOS5"],
        "qv_i":         ["GEOS5"],
        "wv_i":         ["GEOS5"],
        "p_air_i":      ["GEOS5"],
        "p_air_0_i":    ["GEOS5"],

        # Temporal constants
        "lw_offset":    ["STATICS"],
        "lw_slope":     ["STATICS"],
        "r0_bare":      ["STATICS"],
        "r0_full":      ["STATICS"],
        "rn_offset":    ["STATICS"],
        "rn_slope":     ["STATICS"],
        "t_amp_year":   ["STATICS"],
        "t_opt":        ["STATICS"],
        "vpd_slope":    ["STATICS"],
        "z_oro":        ["STATICS"],

        # Level name
        "level_name": "sideloading",
    }

    # Give product folder.
    extra_source_locations = {
        ("LS7", "ndvi"): r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/NDVI",
        ("LS7", "lst"): r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LST",
        ("LS7", "r0"): r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/ALBEDO",
    }

    main(project_folder, startdate, enddate, latlim, lonlim, level = level,
        extra_source_locations = extra_source_locations)