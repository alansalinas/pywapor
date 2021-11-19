import xarray as xr
import numpy as np
from datetime import datetime as dat
import pywapor.collect as c
import pywapor.general as g
import os
import pandas as pd
from datetime import time as datt
import tqdm
import warnings

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1"):

    print("\n#################")
    print("## PRE_SE_ROOT ##")
    print("#################")

    # Disable Warnings
    warnings.filterwarnings('ignore')

    sdate = dat.strptime(startdate, "%Y-%m-%d").date()
    edate = dat.strptime(enddate, "%Y-%m-%d").date()

    #### Check if selected sources and dates are valid
    levels = g.variables.get_source_level_selections()

    if isinstance(level, dict):
        level_name = "custom"
        if "name" in level.keys():
            level_name = level.pop("name")
        levels[level_name] = level
        level = level_name

    source_selection = levels[level]
    succes = g.tests.check_source_selection(source_selection, sdate, edate)[1]
    assert succes, "invalid source_selection"

    raw_folder = os.path.join(project_folder, "RAW")
    temp_folder = os.path.join(project_folder, "temporary")

    dl_args = (raw_folder, latlim, lonlim, startdate, enddate)

    #### NDVI ####
    print(f"\n#### NDVI ####")
    raw_ndvi_files = list()
    # Order is important! PROBV gets priority over MOD13, and MOD13 over MYD13.
    if "PROBAV" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.PROBAV.PROBAV_S5(*dl_args)[0])
    if "MOD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MOD13.NDVI(*dl_args, remove_hdf = 0))
    if "MYD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MYD13.NDVI(*dl_args))

    cmeta = {
        "composite_type": False,
        "temporal_interp": False,
        "spatial_interp": "nearest",
        "var_name": "NDVI",
        "var_unit": "-",
    }
    ds = g.compositer.main(cmeta, raw_ndvi_files, None, temp_folder, None, lean_output = False)
    ds_ndvi = ds.rename({"band_data": "ndvi"})

    #### LST ####
    print(f"\n#### LST ####")
    raw_lst_files = list()
    if "MOD11" in source_selection["LST"]:
        raw_lst_files.append(c.MOD11.LST(*dl_args))
    if "MYD11" in source_selection["LST"]:
        raw_lst_files.append(c.MYD11.LST(*dl_args))
    ds_lst = combine_lst(raw_lst_files)

    #### DEM ####
    print(f"\n#### DEM ####")
    if "SRTM" in source_selection["DEM"]:
        raw_dem_file = c.SRTM.DEM(*dl_args[:3])

    #### METEO ####
    print(f"\n#### METEO ####")
    if "MERRA2" in source_selection["METEO"]:
        freq = "H"
        periods, period_times = calc_periods(ds_lst.time.values, freq)
        ds_lst["periods"] = xr.DataArray([x[2] for x in periods], {"time": ds_lst.time})
        meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
        waitbar = tqdm.tqdm(total = len(periods) * len(meteo_vars), delay = 10, initial = 1)
        for sd, ed, period in periods:
            c.MERRA.hourly_MERRA2(*dl_args[:3], sd, ed, meteo_vars, [int(period)])
    elif "GEOS5" in source_selection["METEO"]:
        freq = "3H"
        periods, period_times = calc_periods(ds_lst.time.values, freq)
        ds_lst["periods"] = xr.DataArray([x[2] for x in periods], {"time": ds_lst.time})
        meteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']
        waitbar = tqdm.tqdm(total = len(periods) * len(meteo_vars), delay = 10, initial = 1)
        for sd, ed, period in periods:
            c.GEOS.three_hourly(*dl_args[:3], sd, ed, meteo_vars, [int(period)], Waitbar = waitbar)

    #### INST. METEO ####
    raw_inst_meteo_paths = g.variables.get_raw_meteo_paths("inst")
    ds_temperature = None
    ds_meteo = None
    for var_name, folder in raw_inst_meteo_paths[source_selection["METEO"][0]].items():
        print(f"## {var_name} ##")
        for date, period in zip(ds_lst.time.values, ds_lst.periods.values):
            hour = int((period - 1) * {"3H": 3, "H": 1}[freq])
            hour_str = str(hour).zfill(2)
            date_str = pd.Timestamp(date).strftime("%Y.%m.%d")
            raw_file = os.path.join(*folder).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
            if var_name == "tair_inst":
                raw_file = g.lapse_rate.lapse_rate_temperature(raw_file, raw_dem_file)
                meteo_datetime = dat.combine(pd.Timestamp(date).date(), period_times[period])
                da_var = xr.open_dataset(raw_file).squeeze("band").rename({"band_data": var_name, "x": "lon", "y": "lat"})[var_name]
                da_var = da_var.assign_coords({"time": meteo_datetime}).expand_dims("time", axis = 0)
                ds_temp = xr.Dataset({var_name: da_var})
                if isinstance(ds_temperature, type(None)):
                    ds_temperature = ds_temp
                elif np.datetime64(meteo_datetime) in ds_temperature.time:
                    ds_temperature = xr.merge([ds_temperature, ds_temp])
                else:
                    ds_temperature = xr.concat([ds_temperature, ds_temp], dim = "time")
            else:
                meteo_datetime = dat.combine(pd.Timestamp(date).date(), period_times[period])
                da_var = xr.open_dataset(raw_file).squeeze("band").rename({"band_data": var_name, "x": "lon", "y": "lat"})[var_name]
                da_var = da_var.assign_coords({"time": meteo_datetime}).expand_dims("time", axis = 0)
                ds_temp = xr.Dataset({var_name: da_var})
                if isinstance(ds_meteo, type(None)):
                    ds_meteo = ds_temp
                elif np.datetime64(meteo_datetime) in ds_meteo.time:
                    ds_meteo = xr.merge([ds_meteo, ds_temp])
                else:
                    ds_meteo = xr.concat([ds_meteo, ds_temp], dim = "time")

    print("\n#################")
    print("## PRE_SE_ROOT ##")
    print("##### DONE ######\n")

    return ds_lst, ds_meteo, ds_ndvi, ds_temperature

def calc_periods(times, freq):
    periods = list()

    for now_time in times:

        date = pd.Timestamp(now_time).date()
        starttime = dat(date.year, date.month, date.day, 0, 0)
        endtime = dat(date.year, date.month, date.day, 23, 59)

        offsets = {"3H": 90, "H": 30}

        DateTime = pd.date_range(starttime, endtime, freq=freq) + pd.offsets.Minute(offsets[freq])
        Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - now_time))
        period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1

        periods.append((starttime, endtime, period))

    period_times = {i+1: datt(time.hour, time.minute) for i, time in enumerate(DateTime)}

    return periods, period_times

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