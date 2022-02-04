# -*- coding: utf-8 -*-
"""
"""
import os
import glob
import warnings
import shutil
import pywapor
import xarray as xr
import pandas as pd
import numpy as np
import pywapor.general as g
import pywapor.general.processing_functions as pf
import pywapor.general.pre_defaults as defaults
from datetime import datetime as dat
from pywapor.collect.downloader import collect_sources
from pywapor.general.logger import log, adjust_logger
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.compositer import calculate_ds
import copy

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1", 
        diagnostics = None, composite_length = "DEKAD", extra_composite_enhancements = {},
        extra_source_enhancements = {}, extra_source_locations = None, 
        se_root_version = "v2", process_vars = "all"):

    # Create project folder.
    if not os.path.exists(project_folder):
        os.makedirs(project_folder)

    # Initiate logger.
    log_write = True
    log_level = "INFO"
    adjust_logger(log_write, project_folder, log_level)

    log.info("> PRE_ET_LOOK").add()

    # Disable Warnings
    warnings.filterwarnings('ignore')

    # Calculate relevant dates.
    sdate = dat.strptime(startdate, "%Y-%m-%d").date()
    edate = dat.strptime(enddate, "%Y-%m-%d").date()
    epochs_info = create_dates(sdate, edate, composite_length)

    # Load required variable sources.
    if isinstance(level, str):
        levels = g.variables.get_source_level_selections()
        source_selection = levels[level]
        level_name = level
    elif isinstance(level, dict):
        source_selection = copy.copy(level)
        if "level_name" in source_selection.keys():
            level_name = source_selection.pop("level_name")
        else:
            level_name = "custom"

    # Define folders.
    level_folder, temp_folder, raw_folder = get_folders(project_folder, level_name)
    if isinstance(diagnostics, dict):
        diagnostics["folder"] = os.path.join(level_folder, "graphs")

    # Define download arguments.
    dl_args = {
                "Dir": raw_folder, 
                "latlim": latlim,
                "lonlim": lonlim,
                "Startdate": startdate, 
                "Enddate": enddate,
                }

    # List of variables to process. Resampling resolution is taken from 
    # first item in the list.
    if process_vars == "all":
        all_vars = [
                    "ndvi",
                    "p_24",
                    "se_root",
                    "r0", 
                    "z",
                    "lulc",
                    "ra_24",
                    't_air_24', 
                    't_air_min_24', 't_air_max_24', 
                    'u2m_24', 'v2m_24', 'p_air_0_24', 'qv_24',
                    'lw_offset', 'lw_slope', 'r0_bare', 'r0_full', 'rn_offset', 
                    'rn_slope', 't_amp_year', 't_opt', 'vpd_slope', 'z_oro',
                    # # 'land_mask', 'rs_min', 'z_obst_max', # these are generated from lulc
                    ]
    else:
        all_vars = process_vars

    datasets = list()

    # Loop over the variables, resampling, enhancing and compositing each one.
    for var in all_vars:

        log.info(f"# {var}")

        # Run pre_se_root and se_root.
        if var == "se_root":
            ds_in = pywapor.pre_se_root.main(project_folder, startdate, enddate, latlim, 
                                            lonlim, level = level,
                                            extra_source_locations = extra_source_locations)[0]
            raw_files = [pywapor.se_root.main(ds_in, se_root_version = se_root_version, 
                                                export_to_tif = True, export_vars = "default")["se_root"]]
        # Download RAW-data.
        else:
            raw_files = collect_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)

        # Define resampling parameters.
        if 'example_info' not in vars():
            example_info = pf.select_template(raw_files)
            # example_info[1].to_netcdf(r"/Users/hmcoerver/pywapor_notebooks/example_ds.nc")

        # Create composites.
        ds = unraw_all(var, raw_files, epochs_info, temp_folder, 
                        example_info[1], diagnostics = diagnostics, 
                        extra_source_enhancements = extra_source_enhancements)

        # Store variable composites in a list.
        datasets.append(ds)

    # Merge all the variable composites into one xr.Dataset.
    ds = xr.merge(datasets, combine_attrs = "drop_conflicts")
    ds.attrs = {}

    log.info("> Composite enhancers.").add()

    # Rename variables.
    ds = ds.rename_vars({x: x.replace("_composite", "") for x in list(ds.variables)})

    # Define composite enhancements.
    composite_enhancements = defaults.composite_enhancements_defaults()
    composite_enhancements = {**composite_enhancements, **extra_composite_enhancements}

    # Apply composite enhancements.
    for variable, enhancers in composite_enhancements.items():
        if variable not in list(ds.keys()):
            continue
        log.info(f"# {variable}")
        for enhancer in enhancers:
            ds, label = apply_enhancer(ds, variable, enhancer)
            ds, _ = calculate_ds(ds, label = label)

    log.sub().info("< Composite enhancers.")

    ds = g.variables.fill_attrs(ds)

    ds.attrs["geotransform"] = example_info[2][0]
    ds.attrs["projection"] = example_info[2][1]
    ds.attrs["pixel_size"] = example_info[3]
    ds.attrs["example_file"] = example_info[0]

    if "spatial_ref" in list(ds.coords):
        ds = ds.drop_vars("spatial_ref")

    # Save pre_et_look-output (i.e. et_look-input).
    all_vars = [var for var in list(ds.variables) if "lon" in ds[var].coords 
                                                and "lat" in ds[var].coords]
    out_fh = os.path.join(level_folder, "et_look_input.nc")
    ds, out_fh = calculate_ds(ds, out_fh, label = "--> Saving results.", 
                                encoding = {k: {"dtype": "float32"} for k in all_vars}
                                )

    # Remove temporary files.
    for fh in glob.glob(os.path.join(temp_folder, "*.nc")):
        os.remove(fh)
    if len(os.listdir(temp_folder)) == 0:
        shutil.rmtree(temp_folder)

    # Reset working directory.
    os.chdir(project_folder)

    log.sub().info("< PRE_ET_LOOK")

    return ds, out_fh

def get_folders(project_folder, level_name = None):
    """Define some folders based on a root-folder.

    Parameters
    ----------
    project_folder : str
        Path to root-folder.
    level_name : str, optional
        Name of the source selection level, by default None.

    Returns
    -------
    str
        Path to level-folder.
    str
        Path to temporary-data-folder.
    str
        Path to raw-folder.
    """
    if not isinstance(level_name, type(None)):
        level_folder = os.path.join(project_folder, level_name)
    else:
        level_folder = None
    temp_folder = os.path.join(project_folder, "temporary")
    raw_folder = os.path.join(project_folder, "RAW")
    return level_folder, temp_folder, raw_folder

def unraw_all(var, raw_files, epochs_info, temp_folder, 
                example_ds, diagnostics = None, extra_source_enhancements = {}):

    cmeta = defaults.composite_defaults()[var]

    ds = g.compositer.main(cmeta, raw_files, epochs_info, temp_folder, example_ds,
                            extra_source_enhancements = extra_source_enhancements)
    if isinstance(diagnostics, dict):
        _ = g.compositer.main(cmeta, raw_files, epochs_info, None, example_ds, 
                                    lean_output = False, diagnostics = diagnostics, 
                                    extra_source_enhancements = extra_source_enhancements)

    return ds

def create_dates(sdate, edate, period_length):
    """Define composite lenghts

    Parameters
    ----------
    sdate : str or datetime.datetime
        Date from which to start.
    edate : str or datetime.datetime
        Enddate.
    period_length : int or str
        Amount of days between the returned dates, can also be "DEKAD" in which
        case it will split each month in roughly 10-day segments.

    Returns
    -------
    np.ndarray
        [description]
    """
    if isinstance(sdate, str):
        sdate = dat.strptime(sdate, "%Y-%m-%d")
    if isinstance(edate, str):
        edate = dat.strptime(edate, "%Y-%m-%d")
    if period_length == "DEKAD":
        dekad1 = pd.date_range(sdate - pd.to_timedelta(sdate.day - 1, unit='d'), edate, freq = "MS")
        dekad2 = dekad1 + pd.Timedelta(10, unit='D')
        dekad3 = dekad1 + pd.Timedelta(20, unit='D')
        dates = np.sort(np.array([dekad1, dekad2, dekad3]).flatten())
        no_periods = dates.size - 1
    else:
        days = (edate - sdate).days
        no_periods = int(days // period_length + 1)
        dates = pd.date_range(sdate, periods = no_periods + 1 , freq = f"{period_length}D")
    periods_start = dates[:-1]
    mask = np.all([periods_start < np.datetime64(edate), 
                   periods_start >= np.datetime64(sdate)], axis = 0)
    epochs = np.arange(0, no_periods, 1)[mask]
    periods_end = dates[1:][mask]
    periods_start = periods_start[mask]
    return epochs, periods_start, periods_end
         
if __name__ == "__main__":

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"
    # startdate = "2021-07-06"
    # enddate = "2021-07-07"
    composite_length = "DEKAD"
    # composite_length = 1
    # level = "level_1"
    extra_sources = None
    extra_source_locations = None
    se_root_version = "v2"

    diagnostics = { # label          # lat      # lon
                    "water":	    (29.44977,	30.58215),
                    "desert":	    (29.12343,	30.51222),
                    "agriculture":	(29.32301,	30.77599),
                    "urban":	    (29.30962,	30.84109),
                    }
    # diagnostics = None

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

    ds_in, fh_in = main(project_folder, startdate, enddate, latlim, lonlim, level = level, 
        diagnostics = diagnostics, composite_length = composite_length,
        extra_source_locations = extra_source_locations, se_root_version = se_root_version)

    # fh_in = r"/Users/hmcoerver/pywapor_notebooks/level_1/et_look_input___.nc"

    # out = pywapor.et_look.main(fh_in, export_vars = "all")

    # param = "t_air_24"
    # sources = ["GEOS5"]

    # dl_args = {
    #             "Dir": os.path.join(project_folder, "RAW"), 
    #             "latlim": latlim, 
    #             "lonlim": lonlim, 
    #             "Startdate": startdate, 
    #             "Enddate": enddate,
    #             }

    # files = collect_sources(param, sources, dl_args)



