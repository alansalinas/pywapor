# -*- coding: utf-8 -*-
"""
"""
import os
from datetime import datetime as dat
from datetime import time as datt
from osgeo import gdal
import pandas as pd
import numpy as np
import requests
import pywapor
from pywapor.et_look import get_geoinfo
import pywapor.general.processing_functions as PF
import pywapor.collect as c
import pywapor.general as g
from pathlib import Path
import json
import xarray as xr

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1", 
        diagnostics = None):

#%%
    sdate = dat.strptime(startdate, "%Y-%m-%d").date()
    edate = dat.strptime(enddate, "%Y-%m-%d").date()
    epochs_info = create_dates(sdate, edate, "DEKAD")

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
    level_folder = os.path.join(project_folder, level)
    graph_folder = os.path.join(level_folder, "graphs")

    dl_args = (raw_folder, latlim, lonlim, startdate, enddate)

    unraw_file_templates = unraw_filepaths(epochs_info[1], level_folder, "{var}")

    #### NDVI #### 
    raw_ndvi_files = list()
    # Order is important! PROBV gets priority over MOD13, and MOD13 over MYD13.
    if "PROBAV" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.PROBAV.PROBAV_S5(*dl_args)[0])
    if "MOD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MOD13.NDVI(*dl_args))
    if "MYD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MYD13.NDVI(*dl_args))

    example_fh, example_ds, example_geoinfo = select_template(raw_ndvi_files)

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "NDVI",
        "var_unit": "-",
    }
    ds_ndvi = compositer(cmeta, raw_ndvi_files, epochs_info, example_ds).rename({"band_data": "ndvi"})
    unraw_all_new(cmeta["var_name"], ds_ndvi, unraw_file_templates, example_geoinfo)

    if not isinstance(diagnostics, type(None)):
        pywapor.post_et_look.plot_composite(ds_ndvi, diagnostics, graph_folder, "ndvi")

    #### ALBEDO ####
    raw_albedo_files = list()
    # Order is important! PROBV gets priority over MDC43.
    if "PROBAV" in source_selection["ALBEDO"]:
        raw_albedo_files.append(c.PROBAV.PROBAV_S5(*dl_args)[1])
    if "MDC43" in source_selection["ALBEDO"]:
        raw_albedo_files.append(c.MCD43.ALBEDO(*dl_args))

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "ALBEDO",
        "var_unit": "-",
    }
    ds = compositer(cmeta, raw_albedo_files, epochs_info, example_ds)
    unraw_all_new(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

    if not isinstance(diagnostics, type(None)):
        pywapor.post_et_look.plot_composite(ds, diagnostics, graph_folder)

    #### LST ####
    raw_lst_files = list()
    if "MOD11" in source_selection["LST"]:
        raw_lst_files.append(c.MOD11.LST(*dl_args))
    if "MYD11" in source_selection["LST"]:
        raw_lst_files.append(c.MYD11.LST(*dl_args))
    
    ds_lst = combine_lst_new(raw_lst_files)

    #### PRECIPITATION ####
    raw_precip_files = list()
    if "CHIRPS" in source_selection["PRECIPITATION"]:
        raw_precip_files.append(c.CHIRPS.daily(*dl_args))

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "Precipitation",
        "var_unit": "mm/day",
    }
    ds = compositer(cmeta, raw_precip_files, epochs_info, example_ds)
    unraw_all_new(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

    if not isinstance(diagnostics, type(None)):
        pywapor.post_et_look.plot_composite(ds, diagnostics, graph_folder)

    #### DEM ####
    if "SRTM" in source_selection["DEM"]:
        raw_dem_file = c.SRTM.DEM(*dl_args[:3])
    
    dem_file = unraw_filepaths(None, level_folder, "DEM", static = True)[0]
    unraw(raw_dem_file, dem_file, example_fh, 4)  

    #### SLOPE ASPECT ####
    slope_aspect(dem_file, level_folder, example_fh)

    #### LULC ####
    if "GLOBCOVER" in source_selection["LULC"]:
        raw_lulc_file = c.Globcover.Landuse(*dl_args[:3])
        raw_lulc_files = [(year, raw_lulc_file) for year in range(sdate.year, edate.year + 1)]
    elif "WAPOR" in source_selection["LULC"]:
        raw_lulc_files = c.WAPOR.Get_Layer(*dl_args[:3], sdate.strftime("%Y-01-01"), edate.strftime("%Y-12-31"), 'L1_LCC_A')
        years = np.array([dat.strptime(os.path.split(fp)[-1], "L1_LCC_A_WAPOR_YEAR_%Y.%m.%d.tif").year for fp in raw_lulc_files])
        raw_lulc_files = [(year, raw_lulc_files[np.argmin(np.abs(years - year))]) for year in range(sdate.year, edate.year + 1)]        

    lulc_values = g.landcover_converter.get_lulc_values()

    lulc_file_template = unraw_filepaths(None, level_folder, "{var}_{year}", static = True)[0]
    for year, raw_file in raw_lulc_files:
        for key, replace_values in lulc_values[source_selection["LULC"][0]].items():
            print(year, raw_file, key)
            unraw_replace_values(raw_file, lulc_file_template.format(var = key, year = year), replace_values, example_fh)

    #### METEO ####
    if "MERRA2" in source_selection["METEO"]:
        freq = "H"
        periods, period_times = calc_periods_new(ds_lst.time.values, freq)
        ds_lst["periods"] = xr.DataArray([x[2] for x in periods], {"time": ds_lst.time})
        meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
        raw_meteo_files = c.MERRA.daily_MERRA2(*dl_args, meteo_vars)
        raw_meteo_files = {**c.MERRA.daily_MERRA2(*dl_args, ['t2m'], data_type = ["min", "max"]), **raw_meteo_files}
        for sd, ed, period in periods:
            c.MERRA.hourly_MERRA2(*dl_args[:3], sd, ed, meteo_vars, [int(period)])
        meteo_name_convertor = {'t2m-max': "tair_max_24", 't2m-min': "tair_min_24", 't2m': "tair_24", 
                                'u2m': "u2m", 'v2m': "v2m",'q2m': "qv_24",'slp': "Pair_24_0"}
    elif "GEOS5" in source_selection["METEO"]:
        freq = "3H"
        periods, period_times = calc_periods_new(ds_lst.time.values, freq)
        ds_lst["periods"] = xr.DataArray([x[2] for x in periods], {"time": ds_lst.time})
        meteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']
        raw_meteo_files = c.GEOS.daily(*dl_args, meteo_vars)
        meteo_name_convertor = {'t2m': "tair_24", 't2m-min': "tair_min_24", 't2m-max': "tair_max_24",  
                                'u2m': "u2m", 'v2m': "v2m",'slp': "Pair_24_0",'qv2m': "qv_24",}
        for sd, ed, period in periods:
            c.GEOS.three_hourly(*dl_args[:3], sd, ed, meteo_vars, [int(period)])

#%%

    #### INST. METEO ####
    raw_inst_meteo_paths = g.variables.get_raw_meteo_paths("inst")
    ds_temperature = None
    ds_meteo = None
    for var_name, folder in raw_inst_meteo_paths[source_selection["METEO"][0]].items():
        for date, period in zip(ds_lst.time.values, ds_lst.periods.values):
            hour = int((period - 1) * {"3H": 3, "H": 1}[freq])
            hour_str = str(hour).zfill(2)
            date_str = pd.Timestamp(date).strftime("%Y.%m.%d")
            raw_file = os.path.join(*folder).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
            if var_name == "tair_inst":
                raw_file = lapse_rate_temp(raw_file, dem_file)
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

#%%

    #### NON-INST. METEO ####
    for var_name, raw_files in raw_meteo_files.items():

        if "t2m" in var_name:
            raw_files = [lapse_rate_temp(x, dem_file) for x in raw_files if "_K_" in x]

        if var_name not in meteo_name_convertor.keys():
            continue

        units = {"tair_max_24": "°C", "tair_min_24":"°C", "tair_24":"°C", 
          "u2m":"m/s", "v2m":"m/s","qv_24":"kg/kg","Pair_24_0":"kPa"}

        cmeta = {
            "composite_type": "mean",
            "temporal_interp": "linear",
            "spatial_interp": "linear",
            "var_name": meteo_name_convertor[var_name],
            "var_unit": units[meteo_name_convertor[var_name]],
        }
        ds = compositer(cmeta, [raw_files], epochs_info, example_ds)
        unraw_all_new(meteo_name_convertor[var_name], ds, unraw_file_templates, example_geoinfo)

        if not isinstance(diagnostics, type(None)):
            pywapor.post_et_look.plot_composite(ds, diagnostics, graph_folder)

    #### LAT LON ####
    lat_lon(level_folder, example_fh)

    #### SOLAR RADIATION ####
    if "MERRA2" in source_selection["TRANS"]:
        raw_ra24_files = c.MERRA.daily_MERRA2(*dl_args, ['swgnet'])['swgnet']

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "linear",
        "var_name": "ra_24",
        "var_unit": "-",
    }
    ds = compositer(cmeta, [raw_ra24_files], epochs_info, example_ds)
    unraw_all_new(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

    if not isinstance(diagnostics, type(None)):
        pywapor.post_et_look.plot_composite(ds, diagnostics, graph_folder)

    #### TEMP. AMPLITUDE ####
    raw_temp_ampl_file = os.path.join(raw_folder, "GLDAS", "Temp_Amplitudes_global.tif")
    download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", raw_temp_ampl_file)
    temp_ampl_file_template = unraw_filepaths(None, level_folder, "Tair_amp_{year}", static = True)[0]
    raw_temp_ampl_files = [(year, raw_temp_ampl_file) for year in range(sdate.year, edate.year + 1)]
    for year, raw_file in raw_temp_ampl_files:
        unraw(raw_file, temp_ampl_file_template.format(year = year), example_fh, 6)

    #### SE_ROOT ####
    raw_se_root_files = calc_se_root_i(project_folder, ds_lst, ds_meteo, ds_ndvi, 
                                        ds_temperature, example_ds, example_geoinfo)

    cmeta = {
        "composite_type": 0.90,
        "temporal_interp": False,
        "spatial_interp": "nearest",
        "var_name": "se_root",
        "var_unit": "-",
    }
    ds = compositer(cmeta, [raw_se_root_files], epochs_info, example_ds)
    unraw_all_new(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

    if not isinstance(diagnostics, type(None)):
        pywapor.post_et_look.plot_composite(ds, diagnostics, graph_folder)

    #### METADATA ####
    metadata = dict()
    metadata["pywapor_version"] = pywapor.__version__
    metadata["created"] = dat.now().strftime("%m/%d/%Y, %H:%M:%S")
    metadata["template_file"] = example_fh
    metadata["geotransform"] = "[{0}, {1}, {2}, {3}, {4}, {5}]".format(*PF.get_geoinfo(example_fh)[0])
    metadata["resolution"] = "[{0}, {1}]".format(*PF.get_geoinfo(example_fh)[2:])
    metadata["inputs"] = dict()
    metadata["inputs"]["latlim"] = "[{0}, {1}]".format(*latlim)
    metadata["inputs"]["lonlim"] = "[{0}, {1}]".format(*lonlim)
    metadata["inputs"]["startdate"] = startdate
    metadata["inputs"]["enddate"] = enddate
    metadata["inputs"]["source_selection_name"] = level
    metadata["inputs"]["project_folder"] = project_folder
    metadata["sources"] = source_selection
    json_file = os.path.join(level_folder, f"metadata_{level}.json")
    with open(json_file, 'w+') as f:
        json.dump(metadata, f, indent = 4 )

    os.chdir(project_folder)

def unraw_all(variable, unraw_file_templates, raw_files, template_file, method):
    unraw_files = [unraw_file.format(var = variable) for unraw_file in unraw_file_templates]
    if len(raw_files) != len(unraw_files):
        matched_raw_files = match_to_nearest(raw_files, unraw_files)
    else:
        matched_raw_files = raw_files
    for raw_file, unraw_file in zip(matched_raw_files, unraw_files):
        unraw(raw_file, unraw_file, template_file, method)
    return unraw_files

def unraw_all_new(var_name, ds, unraw_file_templates, example_geoinfo):
    assert ds.epoch.size == len(unraw_file_templates)
    for i, fh in enumerate(unraw_file_templates):
        fh_date = pd.Timestamp(dat.strptime(fh.split("_")[-1], "%Y%m%d.tif"))
        assert fh_date == pd.Timestamp(ds.epoch_starts.isel(epoch = i).values)
        array = np.copy(ds.composite.isel(epoch = i).values)
        PF.Save_as_tiff(fh.format(var = var_name), 
                        array, example_geoinfo[0], example_geoinfo[1])

def match_to_nearest(raw_files, files):
    raw_dates = [dat.strptime(os.path.split(fp)[-1].split("_")[-1], "%Y.%m.%d.tif") for fp in raw_files]
    dates = [dat.strptime(os.path.split(fp)[-1].split("_")[-1], "%Y%m%d.tif") for fp in files]
    find_idx = lambda d: np.argmin([np.abs((raw_date - d).days) for raw_date in raw_dates])
    matched_raw_ndvi_files = [raw_files[find_idx(d)] for d in dates]
    return matched_raw_ndvi_files

def select_template(fhs):
    fhs = [val for sublist in fhs for val in sublist]

    sizes = [gdal.Open(fh).RasterXSize * gdal.Open(fh).RasterYSize for fh in fhs]
    idx = np.argmax(sizes)

    example_fh = fhs[idx]
    example_ds = xr.open_dataset(example_fh).isel(band = 0).drop_vars(["band", "spatial_ref"]).rename({"x": "lon", "y": "lat"})
    example_geoinfo = PF.get_geoinfo(example_fh)

    return example_fh, example_ds, example_geoinfo
    
def calc_ra24_flat(lat_file, doy):

    ## latitude
    lat = PF.open_as_array(lat_file)

    Gsc = 1367        # Solar constant (W / m2)
    deg2rad = np.pi / 180.0
    # Computation of Hour Angle (HRA = w)
    B = 360./365 * (doy-81)           # (degrees)
    # Computation of cos(theta), where theta is the solar incidence angle
    # relative to the normal to the land surface
    delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)
    
    phi = lat * deg2rad                                     # latitude of the pixel (radians)
    
    dr = 1 + 0.033 * np.cos(doy*2*np.pi/365)
    
    # Daily 24 hr radiation - For flat terrain only !
    ws_angle = np.arccos(-np.tan(phi)*np.tan(delta))   # Sunset hour angle ws   
    
    # Extraterrestrial daily radiation, Ra (W/m2):
    ra24_flat = (Gsc/np.pi * dr * (ws_angle * np.sin(phi) * np.sin(delta) +
                    np.cos(phi) * np.cos(delta) * np.sin(ws_angle)))

    # decl = solar_radiation.declination(doy)
    # iesd = pywapor.et_look_v2.solar_radiation.inverse_earth_sun_distance(doy)
    # ws = pywapor.et_look_v2.solar_radiation.sunset_hour_angle(lat, decl)
    # ra24_flat = pywapor.et_look_v2.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, lat, ws)
    return ra24_flat

def slope_aspect(dem_file, project_folder, template_file):

    slope_file = unraw_filepaths(None, project_folder, "Slope", static = True)[0]
    aspect_file = unraw_filepaths(None, project_folder, "Aspect", static = True)[0]

    if not os.path.exists(slope_file) or not os.path.exists(aspect_file):

        dem = PF.open_as_array(dem_file)

        # constants
        geo_ex, proj_ex, size_x_ex, size_y_ex = PF.get_geoinfo(template_file)
        dlat, dlon = PF.calc_dlat_dlon(geo_ex, size_x_ex, size_y_ex)            

        pixel_spacing = (np.nanmean(dlon) + np.nanmean(dlat)) / 2
        rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree

        # The slope output raster map contains slope values, stated in degrees of inclination from the horizontal
        # Calculate slope
        x, y = np.gradient(dem, pixel_spacing, pixel_spacing)
        hypotenuse_array = np.hypot(x,y)
        slope = np.arctan(hypotenuse_array) * rad2deg

        # calculate aspect
        aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
        aspect = 180 + aspect

        # Save as tiff files
        PF.Save_as_tiff(slope_file, slope, geo_ex, proj_ex)
        PF.Save_as_tiff(aspect_file, aspect, geo_ex, proj_ex)

    return slope_file, aspect_file

def lat_lon(project_folder, template_file):

    lat_file = unraw_filepaths(None, project_folder, "Lat", static = True)[0]
    lon_file = unraw_filepaths(None, project_folder, "Lon", static = True)[0]

    geo_ex, proj_ex, size_x_ex, size_y_ex = PF.get_geoinfo(template_file)

    lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)
    lat_deg = np.array([geo_ex[3] + np.arange(0,size_y_ex) * geo_ex[5]]*size_x_ex).transpose()

    if not os.path.exists(lon_file):
        PF.Save_as_tiff(lon_file, lon_deg, geo_ex, proj_ex)
    if not os.path.exists(lat_file):
        PF.Save_as_tiff(lat_file, lat_deg, geo_ex, proj_ex)

    return lat_file, lon_file

def calc_periods(time_files, freq):

    geo_ex, proj_ex, size_x_ex, size_y_ex = PF.get_geoinfo(time_files[0])
    periods = list()

    for time_file in time_files:

        array = PF.open_as_array(time_file)
        dtime = np.nanmean(array)
        if np.isnan(dtime):
            dtime = 12

        lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)

        offset_GTM = int(round(lon_deg[int(lon_deg.shape[0]/2),int(lon_deg.shape[1]/2)] * 24 / 360))

        date = dat.strptime(os.path.split(time_file)[-1], "Time_%Y%m%d.tif")
        starttime = dat(date.year, date.month, date.day, 0, 0)
        endtime = dat(date.year, date.month, date.day, 23, 59)
        nowtime = dat(date.year, date.month, date.day, int(np.floor(dtime)), int((dtime - np.floor(dtime))*60)) - pd.DateOffset(hours = offset_GTM) 

        offsets = {"3H": 90, "H": 30}

        DateTime = pd.date_range(starttime, endtime, freq=freq) + pd.offsets.Minute(offsets[freq])
        Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - nowtime))
        period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1

        periods.append((starttime, endtime, period))

    return periods

def unraw_filepaths(periods_start, project_folder, var, static = False):
    filepaths = list()
    if not static:
        for date in periods_start:
            date_str = pd.Timestamp(date).strftime("%Y%m%d")
            fp = os.path.join(project_folder,
                            date_str, f"{var}_{date_str}.tif")
            filepaths.append(fp)
    else:
        fp = os.path.join(project_folder,
                        "static", f"{var}.tif")
        filepaths.append(fp)   
    return filepaths

def create_dates(sdate, edate, period_length):
    if period_length == "DEKAD":
        dekad1 = pd.date_range(sdate, edate, freq = "MS")
        dekad2 = dekad1 + pd.Timedelta(10, unit='D')
        dekad3 = dekad1 + pd.Timedelta(20, unit='D')
        dates = np.sort(np.array([dekad1, dekad2, dekad3]).flatten())
        no_periods = dates.size - 1
    else:
        days = (edate - sdate).days
        no_periods = int(days // period_length + 1)
        dates = pd.date_range(sdate, periods = no_periods + 1 , freq = f"{period_length}D")
    
    periods_start = dates[:-1]
    epochs = np.arange(0, no_periods, 1)[periods_start < np.datetime64(edate)]
    periods_end = np.array([x - np.timedelta64(1, "D") for x in dates[1:]])
    periods_end = periods_end[periods_start < np.datetime64(edate)]
    periods_start = periods_start[periods_start < np.datetime64(edate)]

    return epochs, periods_start, periods_end

def unraw(raw_file, unraw_file, template_file, method):
    if not os.path.exists(unraw_file) and os.path.exists(raw_file):
        geo_ex, proj_ex = PF.get_geoinfo(template_file)[0:2]
        array = PF.reproj_file(raw_file, template_file, method)
        PF.Save_as_tiff(unraw_file, array, geo_ex, proj_ex)

def unraw_replace_values(raw_file, unraw_file, replace_values, template_file):

    if not os.path.exists(unraw_file) and os.path.exists(raw_file):

        array = PF.reproj_file(raw_file, template_file, 1)

        replaced_array = np.ones_like(array) * np.nan
        for key, value in replace_values.items():
            replaced_array[array == key] = value

        geo_ex, proj_ex = PF.get_geoinfo(template_file)[0:2]
        PF.Save_as_tiff(unraw_file, replaced_array, geo_ex, proj_ex)

def lapse_rate_temp(tair_file, dem_file):

    # import matplotlib.pyplot as plt

    # tair_file = r"/Volumes/Data/pre_et_look_NEW/RAW/MERRA/Air_Temperature/daily_MERRA2/t2m_MERRA_K_daily_2019.07.06.tif"
    # tair_file = r"/Volumes/Data/pre_et_look_NEW/RAW/GEOS/Air_Temperature/daily/t2m_GEOS_K_daily_2019.07.06.tif"

    # def _plot_array(array):
    #     plt.clf()
    #     plt.imshow(array)
    #     plt.title(array.shape)
    #     plt.colorbar()
    #     plt.show()

    ## 1
    ds_t_down = PF.reproject_dataset_example(tair_file, dem_file, 2)
    tempe = PF.open_as_array(ds_t_down)
    # _plot_array(tempe - 273.15)

    ## 2
    dem_down = PF.open_as_array(dem_file)
    # _plot_array(dem_down)

    ## 3
    ds_dem_up = PF.reproject_dataset_example(dem_file, tair_file, 4)
    
    dem_up = PF.open_as_array(ds_dem_up)
    dem_up[np.isnan(dem_up)] = 0.
    ds_dem_up = PF.Save_as_MEM(dem_up, ds_dem_up.GetGeoTransform(), PF.Get_epsg(ds_dem_up))
    
    ds_dem_up_down = PF.reproject_dataset_example(ds_dem_up, dem_file, 2)
    dem_up_ave = PF.open_as_array(ds_dem_up_down)
    # _plot_array(dem_up_ave)

    ## Correct wrong values
    dem_down[dem_down <= 0] = 0
    dem_up_ave[dem_up_ave <= 0] = 0

    tdown = pywapor.et_look_v2.meteo.disaggregate_air_temperature(tempe, dem_down, dem_up_ave, lapse = pywapor.et_look_v2.constants.lapse)
    # _plot_array(tdown)

    fh = tair_file.replace("_K_", "_C_")
    geo, projection = get_geoinfo(dem_file)[1:]

    PF.Save_as_tiff(fh, tdown, geo, projection)

    # test_tair = r"/Volumes/Data/pre_et_look_ORIGINAL/ETLook_input_MODIS/20190706/tair_24_20190706.tif"
    # tempe_test = PF.open_as_array(test_tair)
    # _plot_array(tempe_test - tdown)

    return fh

def combine_lst(raw_files):

    new_dict = PF.combine_dicts(raw_files)

    lst_files = list()
    time_files = list()

    for date, scenes in new_dict.items():

        date_str = date.strftime("%Y.%m.%d")

        for i, scene in enumerate(scenes):

            lst = PF.open_as_array(scene[0])
            time = PF.open_as_array(scene[1])
            vza = np.abs(PF.open_as_array(scene[2]))

            if i == 0:
                lsts = np.array([lst])
                vzas = np.array([vza])
                times = np.array([time])
            else:
                lsts = np.concatenate((lsts, np.array([lst])))
                vzas = np.concatenate((vzas, np.array([vza])))
                times = np.concatenate((times, np.array([time])))

        idxs = np.argmin(vzas, axis = 0)

        lst = PF.apply_mask(lsts, idxs, axis = 0)
        time = PF.apply_mask(times, idxs, axis = 0)

        template = scene[0]
        geo_ex, proj_ex = PF.get_geoinfo(template)[0:2]

        out_folder = os.path.join(str(Path(template).parent.parent), "LST") 
        out_file_lst = os.path.join(out_folder, f"LST_MCD11A1_K_daily_{date_str}.tif")
        out_file_time = os.path.join(out_folder, f"Time_MCD11A1_hour_daily_{date_str}.tif")

        lst_files.append(out_file_lst)
        time_files.append(out_file_time)

        PF.Save_as_tiff(out_file_lst, lst, geo_ex, proj_ex)
        PF.Save_as_tiff(out_file_time, time, geo_ex, proj_ex)

    return lst_files, time_files

def download_file_from_google_drive(id, destination):
    
    folder = os.path.split(destination)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    URL = "https://docs.google.com/uc?export=download"

    session = requests.Session()

    response = session.get(URL, params = { 'id' : id }, stream = True)
    token = get_confirm_token(response)

    if token:
        params = { 'id' : id, 'confirm' : token }
        response = session.get(URL, params = params, stream = True)

    save_response_content(response, destination)            

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)

def compositer(cmeta, dbs, epochs_info, example_ds = None):
    """
    composite_type = "max", "mean", "min"
    temporal_interp = False, "linear", "nearest", "zero", "slinear", "quadratic", "cubic"
    spatial_interp = "nearest", "linear"
    """

    composite_type = cmeta["composite_type"]
    temporal_interp = cmeta["temporal_interp"]
    spatial_interp = cmeta["spatial_interp"]

    styling = { 1.0: ("*", "r", 1.0, "MOD13Q1"),
                2.0: ("o", "g", 1.0, "MYD13Q1"),
                3.0: ("v", "b", 1.0, "PROBAV"),
                4.0: ("s", "y", 1.0, "MCD43A3"),
                5.0: ("*", "purple", 1.0, "CHIRPS.v2.0"),
                6.0: ("p", "darkblue", 1.0, "GEOS"),
                7.0: ("h", "gray", 1.0, "MERRA"),
                999.0: ("P", "orange", 1.0, "-"),

                0.0: (".", "k", 0.7, "Interp.")}

    epochs, epoch_starts, epoch_ends = epochs_info

    def preprocess_func(ds):

        date_string = ds.encoding["source"].split("_")[-1]
        if len(date_string.split(".")) > 4:
            date = dat.strptime(date_string, '%Y.%m.%d.%H.%M.tif')
        else:
            date = dat.strptime(date_string, '%Y.%m.%d.tif')

        source = ds.encoding["source"].split("_")[-4]

        ds = ds.drop_vars("spatial_ref")
        ds = ds.squeeze("band")
        ds = ds.drop_vars("band")
        ds = ds.expand_dims("time", axis = 0)
        ds = ds.assign_coords({"time": [date]})
        ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
        ds = ds.swap_dims({"x": f"lon", "y": f"lat"})

        sources = {v[3]: k for k, v in styling.items()}

        ds["band_data"] = ds.band_data.expand_dims("source")
        ds = ds.assign_coords({"source": [source]})

        ds["sources"] = np.isfinite(ds.band_data) * sources[source]

        return ds

    ds = None
    dss = list()
    dbs_names = list()

    for db in dbs:
        sub_ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", 
                                    preprocess = preprocess_func)
        dbs_names.append(sub_ds.source.values[0])

        if isinstance(example_ds, type(None)):
            example_ds = sub_ds.isel(time = 0, 
                                    source = 0).drop_vars(["source", "time"])
        else:
            sub_ds = sub_ds.interp_like(example_ds, method = spatial_interp)

        dss.append(sub_ds)

    ds = xr.merge(dss)

    for i, db_name in enumerate(dbs_names):
        if i == 0:
            ds_temp = ds.sel(source = db_name)
        else:
            ds_temp = ds_temp.fillna(ds.sel(source = db_name))

    ds = ds_temp

    if temporal_interp:
        ds_resampled = ds.resample({"time": "1D"}).interpolate(temporal_interp).drop("sources")
        ds_resampled["sources"] = ds.sources
        ds_resampled["sources"] = ds_resampled["sources"].fillna(0.0)
        ds = ds_resampled

    da1 = xr.DataArray(epochs, coords={"time":epoch_starts})
    da2 = xr.DataArray(epochs, coords={"time":epoch_ends})
    da11 = da1.reindex_like(ds["time"], method = "ffill")
    da22 = da2.reindex_like(ds["time"], method = "bfill")
    ds["epochs"] = da11.where(np.isfinite(da22), np.nan)

    if composite_type == "max":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).max().rename({"epochs": "epoch"})
    elif composite_type == "mean":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).mean().rename({"epochs": "epoch"})
    elif composite_type == "min":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).min().rename({"epochs": "epoch"})
    elif isinstance(composite_type, float):
        ds["composite"] = ds.band_data.chunk(dict(time=-1)).groupby(ds["epochs"]).quantile(composite_type).rename({"epochs": "epoch"})
    else:
        print("No valid composite_type selected.")

    ds["epoch_starts"] = xr.DataArray(epoch_starts, coords = {"epoch": epochs})
    ds["epoch_ends"] = xr.DataArray(epoch_ends, coords = {"epoch": epochs})

    checklist = {True: dbs_names + ["Interp."], False: dbs_names}[True] # TODO interp data should not be included when not used
    sources_styling = {str(k): v for k, v in styling.items() if v[3] in checklist}
    ds.attrs = {**cmeta, **sources_styling}

    return ds

def combine_lst_new(lst_files):
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

def calc_periods_new(times, freq):
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

def calc_se_root_i(project_folder, ds_lst, ds_meteo, ds_ndvi, 
                    ds_temperature, example_ds, example_geoinfo):

    # spatial interpolation
    ds_lst2 = ds_lst.interp_like(example_ds, method = "linear") # TODO check meteo download and remove extrapolation
    ds_meteo2 = ds_meteo.interp_like(example_ds, method = "linear")
    
    # temporal interpolation
    ds_meteo3 = ds_meteo2.interp(time = ds_lst2.time, method = "nearest", kwargs={"fill_value": "extrapolate"},) # TODO download more meteo data and switch to linear
    ds_ndvi2 = ds_ndvi.interp(time = ds_lst.time, method = "linear", kwargs={"fill_value": "extrapolate"},)
    ds_temperature2 = ds_temperature.interp(time = ds_lst2.time, method = "nearest", kwargs={"fill_value": "extrapolate"},)

    req_vars = ['time', 'angle', 'lst', 'lon', 'lat', 'Pair_inst_0', 'Pair_inst',
                'qv_inst', 'tair_inst', 'wv_inst', 'v2m_inst', 'u2m_inst', 'ndvi']

    renames = {"Pair_inst_0": "p_air_0_i", "Pair_inst": "p_air_i", # TODO Make this renaming unnesecary
                "tair_inst": "t_air_i", "qv_inst": "qv_i", 
                "lon": "lon_deg", "wv_inst": "wv_i", "lat": "lat_deg"}

    ds_se_root = xr.merge([ds_lst2, ds_meteo3, ds_ndvi2, ds_temperature2])
    ds_se_root = ds_se_root.drop_vars([x for x in list(ds_se_root.variables) if x not in req_vars])
    ds_se_root = ds_se_root.rename({k: v for k, v in renames.items() if k in list(ds_se_root.variables)})
    ds_se_root["u_i"] = np.sqrt(ds_se_root["v2m_inst"]**2 + ds_se_root["u2m_inst"]**2)

    # ds_se_root.to_netcdf(r"/Users/hmcoerver/On My Mac/se_root_full_test.nc")
    # ds_se_root = xr.open_dataset(r"/Users/hmcoerver/On My Mac/se_root_full_test.nc")

    se_root_folder = os.path.join(project_folder, "SMC")
    if not os.path.exists(se_root_folder):
        os.mkdir(se_root_folder)

    files = list()

    for t in ds_se_root.time.values: # TODO adjust pywapor.et_look.se_root to remove this forloop

        id = ds_se_root.sel(time = t)
        date = pd.Timestamp(t)

        dt_str = date.strftime("%Y.%m.%d.%H.%M")
        fn = f"SMC_-_-_inst_{dt_str}.tif"
        fh = os.path.join(se_root_folder, fn)

        if not os.path.exists(fh): # TODO also check if fh has correct geoinfo
            se_root_i = pywapor.et_look.se_root(id, None, pywapor.et_look_v2, date)[0]["se_root"]
            PF.Save_as_tiff(fh, se_root_i.values, example_geoinfo[0], example_geoinfo[1])

        files.append(fh)

    return files

if __name__ == "__main__":

    project_folder = r"/Volumes/Data/FAO/WaPOR_vs_pyWaPOR/pyWAPOR_v1"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-10"
    level = "level_1"
    diagnostics = {"POI #1":  (29.6785, 31.1955),
                    "POI #2": (29.2, 30.621)}