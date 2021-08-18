# -*- coding: utf-8 -*-
"""
"""

import os
import shutil
import datetime
from osgeo import gdal
import pandas as pd
import numpy as np
import rasterio as rio
import requests
import pywapor
import pywapor.functions.Processing_Functions as PF
from pywapor.functions.Swets_Filter import swets_filter
from pathlib import Path
import password as passwords

def check_source_selection(source_selection, startdate, enddate):

    valid_sources, valid_dates = pywapor.general.variables.get_source_validations()

    temporal_sources = ["METEO", "NDVI", "ALBEDO", "LST", 
                    "PRECIPITATION", "TRANS"]

    check_keys = np.all([key in valid_sources.keys() for key in source_selection.keys()])
    assert check_keys, "invalid key in source_selection"

    assert len(source_selection["DEM"]) == 1, "only one DEM source can be selected"
    assert len(source_selection["METEO"]) == 1, "only one METEO source can be selected"
    assert len(source_selection["PRECIPITATION"]) == 1, "only one PRECIPITATION source can be selected"
    assert len(source_selection["LULC"]) == 1, "only one LULC source can be selected"
    assert len(source_selection["TRANS"]) == 1, "only one TRANS source can be selected"

    results = dict()
    all_results = list()

    for var, sources in source_selection.items():

        check1 = [source in valid_sources[var] for source in sources]
        if var in temporal_sources:
            check2 = [startdate >= valid_dates[source][0] for source in sources]
            check3 = [enddate <= valid_dates[source][1] for source in sources]
        else:
            check2 = [True]
            check3 = [True]

        results[var] = {source: {"valid_source:": check1[i],
                                 "valid_startdate:": check2[i],
                                 "valid_enddate:": check3[i],
                                } for i, source in enumerate(sources)}

        all_results.append(np.all([check1, check2, check3]))
    
    succes = np.all(all_results)

    return results, succes

def prepare_et_look_input(project_folder, startdate, enddate, latlim, lonlim, level = "level_1"):

    sdate = datetime.datetime.strptime(startdate, "%Y-%m-%d").date()
    edate = datetime.datetime.strptime(enddate, "%Y-%m-%d").date()

    levels = pywapor.general.variables.get_source_level_selections()
    lulc_values = pywapor.general.landcover_converter.get_lulc_values()
                        
    source_selection = levels[level]
    succes = check_source_selection(source_selection, sdate, edate)[1]
    assert succes, "invalid source_selection"

    raw_folder = os.path.join(project_folder, "RAW")

    un_nasa, pw_nasa = passwords.passes["NASA"]
    wapor_token = passwords.passes["WAPOR"][1]

    #### NDVI ####
    if "MOD13" in source_selection["NDVI"]:
        dt = datetime.timedelta(days = 8)
        mod13_files = pywapor.collect.MOD13.NDVI(raw_folder, sdate - dt, edate + dt, latlim, lonlim, un_nasa, pw_nasa)
    if "MYD13" in source_selection["NDVI"]:
        dt = datetime.timedelta(days = 8)
        myd13_files = pywapor.collect.MYD13.NDVI(raw_folder, sdate - dt, edate + dt, latlim, lonlim, un_nasa, pw_nasa)
    
    raw_ndvi_files = mod13_files + myd13_files
    raw_dates = [datetime.datetime.strptime(os.path.split(fp)[-1].split("_")[-1], "%Y.%m.%d.tif") for fp in raw_ndvi_files]

    ndvi_files = unraw_filepaths(startdate, enddate, project_folder, "NDVI")
    ndvi_dates = [datetime.datetime.strptime(os.path.split(fp)[-1], "NDVI_%Y%m%d.tif") for fp in ndvi_files]
    
    if not os.path.exists(os.path.split(ndvi_files[0])[0]):
        os.makedirs(os.path.split(ndvi_files[0])[0])

    find_idx = lambda d: np.argmin([np.abs((raw_date - d).days) for raw_date in raw_dates])
    idxs = [find_idx(d) for d in ndvi_dates]

    for idx, unraw_file in zip(idxs, ndvi_files):
        shutil.copy(raw_ndvi_files[idx], unraw_file)

    template_file = raw_ndvi_files[0]

    #### ALBEDO ####
    if "MDC43" in source_selection["ALBEDO"]:
        raw_albedo_files = pywapor.collect.MCD43.ALBEDO(raw_folder, startdate, enddate, latlim, lonlim, un_nasa, pw_nasa)

    albedo_files = unraw_filepaths(startdate, enddate, project_folder, "ALBEDO")
    for raw_file, unraw_file in zip(raw_albedo_files, albedo_files):
        unraw(raw_file, unraw_file, template_file, 1)
   
    #### LST ####
    if "MOD11" in source_selection["LST"]:
        raw_mod11_files = pywapor.collect.MOD11.LST(raw_folder, startdate, enddate, latlim, lonlim, un_nasa, pw_nasa)
    if "MYD11" in source_selection["LST"]:
        raw_myd11_files = pywapor.collect.MYD11.LST(raw_folder, startdate, enddate, latlim, lonlim, un_nasa, pw_nasa)
    raw_lst_files, raw_time_files = combine_lst(raw_folder, startdate, enddate)[:2] # TODO: make function flexible

    lst_files = unraw_filepaths(startdate, enddate, project_folder, "LST")
    for raw_file, unraw_file in zip(raw_lst_files, lst_files):
        unraw(raw_file, unraw_file, template_file, 2)

    #### TIME ####
    time_files = unraw_filepaths(startdate, enddate, project_folder, "Time")
    for raw_file, unraw_file in zip(raw_time_files, time_files):
        unraw(raw_file, unraw_file, template_file, 1)

    #### PRECIPITATION ####
    if "CHIRPS" in source_selection["PRECIPITATION"]:
        raw_precip_files = pywapor.collect.CHIRPS.daily(raw_folder, startdate, enddate, latlim, lonlim)
    
    precip_files = unraw_filepaths(startdate, enddate, project_folder, "Precipitation")
    for raw_file, unraw_file in zip(raw_precip_files, precip_files):
        unraw(raw_file, unraw_file, template_file, 6)

    #### DEM ####
    if "SRTM" in source_selection["DEM"]:
        raw_dem_file = pywapor.collect.SRTM.DEM(raw_folder, latlim, lonlim)
    dem_file = unraw_filepaths("", "", project_folder, "DEM", static = True)[0]
    unraw(raw_dem_file, dem_file, template_file, 4)

    #### SLOPE ASPECT ####
    slope_aspect(dem_file, project_folder, template_file)

    #### LULC ####
    if "GLOBCOVER" in source_selection["LULC"]:
        raw_lulc_file = pywapor.collect.Globcover.Landuse(raw_folder, latlim, lonlim)
        raw_lulc_files = [(year, raw_lulc_file) for year in range(sdate.year, edate.year + 1)]
    elif "WAPOR" in source_selection["LULC"]:
        raw_lulc_files = pywapor.collect.WAPOR.Get_Layer(raw_folder, sdate.strftime("%Y-01-01"), edate.strftime("%Y-12-31"), latlim, lonlim, 'L1_LCC_A', wapor_token)
        raw_lulc_files = [(datetime.datetime.strptime(os.path.split(fp)[-1], "L1_LCC_A_WAPOR_YEAR_%Y.%m.%d.tif").year, fp) for fp in raw_lulc_files]

    lulc_file_template = unraw_filepaths("", "", project_folder, "{var}_{year}", static = True)[0]
    for year, raw_file in raw_lulc_files:
        for key, replace_values in lulc_values[source_selection["LULC"][0]].items():
            print(year, raw_file, key)
            unraw_replace_values(raw_file, lulc_file_template.format(var = key, year = year), replace_values, template_file)

    #### METEO ####
    if "MERRA2" in source_selection["METEO"]:
        freq = "H"
        periods = calc_periods(time_files, freq)
        meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
        pywapor.collect.MERRA.daily_MERRA2(raw_folder, meteo_vars, startdate, enddate, latlim, lonlim, un_nasa, pw_nasa)
        pywapor.collect.MERRA.daily_MERRA2(raw_folder, ['t2m'], startdate, enddate, latlim, lonlim, un_nasa, pw_nasa, data_type = ["mean", "min", "max"])
        for sd, ed, period in periods:
            pywapor.collect.MERRA.hourly_MERRA2(raw_folder, meteo_vars, sd, ed, latlim, lonlim, un_nasa, pw_nasa, [int(period)])
    elif "GEOS5" in source_selection["METEO"]:
        freq = "3H"
        periods = calc_periods(time_files, freq)
        meteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']
        pywapor.collect.GEOS.daily(raw_folder, meteo_vars, startdate, enddate, latlim, lonlim)
        for sd, ed, period in periods:
            pywapor.collect.GEOS.three_hourly(raw_folder, meteo_vars, sd, ed, latlim, lonlim, [int(period)])

    meteo_files_template = unraw_filepaths(startdate, enddate, project_folder, "{var}")
    
    raw_meteo_paths = pywapor.general.variables.get_raw_meteo_paths()

    for meteo_file_template in meteo_files_template:
        date = datetime.datetime.strptime(os.path.split(meteo_file_template)[-1], "{var}_%Y%m%d.tif")
        date_str = date.strftime("%Y.%m.%d")
        period = [x[2] for x in periods if x[0] == date][0] 
        hour = int((period - 1) * {"3H": 3, "H": 1}[freq])
        hour_str = str(hour).zfill(2)
        for key, path in raw_meteo_paths[source_selection["METEO"][0]].items():
            unraw_file = meteo_file_template.format(var = key)
            if "wind" in key:
                raw_file_u = os.path.join(*path[0]).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                raw_file_v = os.path.join(*path[1]).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                u_wind = reproj_file(raw_file_u, template_file, 6)
                v_wind = reproj_file(raw_file_v, template_file, 6)
                wind = np.sqrt(u_wind**2 + v_wind**2)
                geo_ex, proj_ex = get_geoinfo(template_file)[0:2]
                PF.Save_as_tiff(unraw_file, wind, geo_ex, proj_ex)
            elif "tair" in key:
                raw_file = os.path.join(*path).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                tair = lapse_rate_temp(raw_file, dem_file)
                geo_ex, proj_ex = get_geoinfo(template_file)[0:2]
                PF.Save_as_tiff(unraw_file, tair, geo_ex, proj_ex)
            else:
                raw_file = os.path.join(*path).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                unraw(raw_file, unraw_file, template_file, method = 6)

    #### LAT LON ####
    lat_file = lat_lon(project_folder, template_file)[0]

    #### TRANS ####
    if "MERRA2" in source_selection["TRANS"]:
        raw_trans_files = pywapor.collect.MERRA.daily_MERRA2(raw_folder, ['swgnet'], startdate, enddate, latlim, lonlim, un_nasa, pw_nasa)[0]
    
    trans_files = unraw_filepaths(startdate, enddate, project_folder, "Trans")
    trans_files = [(int(datetime.datetime.strptime(os.path.split(fp)[-1], "Trans_%Y%m%d.tif").strftime("%j")), fp) for fp in trans_files]

    for (doy, unraw_file), raw_file in zip(trans_files, raw_trans_files):
        ra24_flat = calc_ra24_flat(lat_file, doy)
        array = reproj_file(raw_file, template_file, 6) / ra24_flat
        geo_ex, proj_ex = get_geoinfo(template_file)[0:2]
        PF.Save_as_tiff(unraw_file, array, geo_ex, proj_ex)

    #### TEMP. AMPLITUDE ####
    raw_temp_ampl_file = os.path.join(raw_folder, "GLDAS", "Temp_Amplitudes_global.tif")
    download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", raw_temp_ampl_file)
    temp_ampl_file_template = unraw_filepaths("", "", project_folder, "Tair_amp_{year}", static = True)[0]
    raw_temp_ampl_files = [(year, raw_temp_ampl_file) for year in range(sdate.year, edate.year + 1)]
    for year, raw_file in raw_temp_ampl_files:
        unraw(raw_file, temp_ampl_file_template.format(year = year), template_file, 6)
    
def calc_ra24_flat(lat_file, doy):
    ## latitude
    lat = open_as_array(lat_file)

    ## declination
    deg2rad = np.pi / 180.0
    B = 360./365 * (doy - 81)
    decl = np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B)))
    # decl = solar_radiation.declination(doy)

    iesd = pywapor.et_look_v2.solar_radiation.inverse_earth_sun_distance(doy)
    ws = pywapor.et_look_v2.solar_radiation.sunset_hour_angle(lat, decl)
    ra24_flat = pywapor.et_look_v2.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, lat, ws)
    return ra24_flat

def slope_aspect(dem_file, project_folder, template_file):

    slope_file = unraw_filepaths("", "", project_folder, "Slope", static = True)[0]
    aspect_file = unraw_filepaths("", "", project_folder, "Aspect", static = True)[0]

    if not os.path.exists(slope_file) or not os.path.exists(aspect_file):

        dem = open_as_array(dem_file)

        # constants
        geo_ex, proj_ex, size_x_ex, size_y_ex = get_geoinfo(template_file)
        dlat, dlon = calc_dlat_dlon(geo_ex, size_x_ex, size_y_ex)            

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

    lat_file = unraw_filepaths("", "", project_folder, "Lat", static = True)[0]
    lon_file = unraw_filepaths("", "", project_folder, "Lon", static = True)[0]

    geo_ex, proj_ex, size_x_ex, size_y_ex = get_geoinfo(template_file)

    lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)
    lat_deg = np.array([geo_ex[3] + np.arange(0,size_y_ex) * geo_ex[5]]*size_x_ex).transpose()

    if not os.path.exists(lon_file):
        PF.Save_as_tiff(lon_file, lon_deg, geo_ex, proj_ex)
    if not os.path.exists(lat_file):
        PF.Save_as_tiff(lat_file, lat_deg, geo_ex, proj_ex)

    return lat_file, lon_file

def calc_periods(time_files, freq):

    geo_ex, proj_ex, size_x_ex, size_y_ex = get_geoinfo(time_files[0])
    periods = list()

    for time_file in time_files:

        array = open_as_array(time_file)
        dtime = np.nanmean(array)
        if np.isnan(dtime):
            dtime = 12

        lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)

        offset_GTM = int(round(lon_deg[int(lon_deg.shape[0]/2),int(lon_deg.shape[1]/2)] * 24 / 360))

        date = datetime.datetime.strptime(os.path.split(time_file)[-1], "Time_%Y%m%d.tif")
        starttime = datetime.datetime(date.year, date.month, date.day, 0, 0)
        endtime = datetime.datetime(date.year, date.month, date.day, 23, 59)
        nowtime = datetime.datetime(date.year, date.month, date.day, int(np.floor(dtime)), int((dtime - np.floor(dtime))*60)) - pd.DateOffset(hours = offset_GTM) 

        offsets = {"3H": 90, "H": 30}

        DateTime = pd.date_range(starttime, endtime, freq=freq) + pd.offsets.Minute(offsets[freq])
        Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - nowtime))
        period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1

        periods.append((starttime, endtime, period))

    return periods

def unraw_filepaths(startdate, enddate, project_folder, var, static = False):
    filepaths = list()
    if not static:
        dates = pd.date_range(startdate, enddate, freq="D")
        for date in dates:
            date_str = date.strftime("%Y%m%d")
            fp = os.path.join(project_folder, "et_look_input",
                            date_str, f"{var}_{date_str}.tif")
            filepaths.append(fp)
    else:
        fp = os.path.join(project_folder, "et_look_input",
                        "static", f"{var}.tif")
        filepaths.append(fp)   
    return filepaths

def unraw(raw_file, unraw_file, template_file, method):
    if not os.path.exists(unraw_file) and os.path.exists(raw_file):
        geo_ex, proj_ex = get_geoinfo(template_file)[0:2]
        array = reproj_file(raw_file, template_file, method)
        if unraw_file != "":
            PF.Save_as_tiff(unraw_file, array, geo_ex, proj_ex)
        return array

def reproj_file(file, template, method):
    ds = PF.reproject_dataset_example(file, template, method = method)
    array = open_as_array(ds)
    return array

def open_as_array(input):
    if isinstance(input, str):
        ds = gdal.Open(input)
    elif isinstance(input, gdal.Dataset):
        ds = input
    array = ds.GetRasterBand(1).ReadAsArray()
    ndv = ds.GetRasterBand(1).GetNoDataValue()
    array[np.isnan(array)] = ndv
    return array  

def get_geoinfo(template_file):
    ds = gdal.Open(template_file)
    geo_ex = ds.GetGeoTransform()
    proj_ex = ds.GetProjection()
    size_x_ex = ds.RasterXSize
    size_y_ex = ds.RasterYSize
    return (geo_ex, proj_ex, size_x_ex, size_y_ex)

def unraw_replace_values(raw_file, unraw_file, replace_values, template_file):

    if not os.path.exists(unraw_file) and os.path.exists(raw_file):

        array = reproj_file(raw_file, template_file, 1)

        replaced_array = np.ones_like(array) * np.nan
        for key, value in replace_values.items():
            replaced_array[array == key] = value

        geo_ex, proj_ex = get_geoinfo(template_file)[0:2]
        PF.Save_as_tiff(unraw_file, replaced_array, geo_ex, proj_ex)

def lapse_rate_temp(tair_file, dem_file):

    ds_t_down = PF.reproject_dataset_example(tair_file, dem_file, 2)
    ds_dem_up = PF.reproject_dataset_example(dem_file, tair_file, 4)
    ds_dem_up_down = PF.reproject_dataset_example(ds_dem_up, dem_file, 2)

    # Open Arrays
    # T = ds_t_down.GetRasterBand(1).ReadAsArray()
    tempe = open_as_array(ds_t_down)
    # destDEM_down = gdal.Open(dem_file)
    # dem_down = destDEM_down.GetRasterBand(1).ReadAsArray()
    dem_down = open_as_array(dem_file)
    # dem_up_ave = ds_dem_up_down.GetRasterBand(1).ReadAsArray()
    dem_up_ave = open_as_array(ds_dem_up_down)

    # correct wrong values
    dem_down[dem_down <= 0] = 0
    dem_up_ave[dem_up_ave <= 0] = 0

    tdown = pywapor.et_look_v2.meteo.disaggregate_air_temperature(tempe, dem_down, dem_up_ave)

    return tdown

def combine_dicts(dicts):
    new_dict = dict()
    for d in dicts:
        for key, value in d.items():
            if key in new_dict.keys():
                new_dict[key].append(value)
            else:
                new_dict[key] = [value]
    return new_dict

def combine_lst(folders_input_RAW, Startdate, Enddate):

    Dates = pd.date_range(Startdate, Enddate, freq = "D")

    output_folder_end = os.path.join(folders_input_RAW, "MODIS", "LST")

    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)

    lst_files = list()
    time_files = list()
    vza_files = list()

    for Date in Dates:

        LST_file = os.path.join(output_folder_end, "LST_MCD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        Time_file = os.path.join(output_folder_end, "Time_MCD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        VZA_file = os.path.join(output_folder_end, "VZA_MCD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))

        if not (os.path.exists(Time_file) and os.path.exists(LST_file) and os.path.exists(VZA_file)):
            filename_angle_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "Angle_MOD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_angle_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "Angle_MYD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "Time_MOD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "Time_MYD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_lst_mod = os.path.join(folders_input_RAW, "MODIS", "MOD11", "LST_MOD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_lst_myd = os.path.join(folders_input_RAW, "MODIS", "MYD11", "LST_MYD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))

            dest_angle_mod = gdal.Open(filename_angle_mod)
            dest_angle_myd = gdal.Open(filename_angle_myd)
            dest_time_mod = gdal.Open(filename_time_mod)
            dest_time_myd = gdal.Open(filename_time_myd)
            dest_lst_mod = gdal.Open(filename_lst_mod)
            dest_lst_myd = gdal.Open(filename_lst_myd)

            Array_angle_mod = dest_angle_mod.GetRasterBand(1).ReadAsArray()
            Array_angle_myd = dest_angle_myd.GetRasterBand(1).ReadAsArray()
            Array_time_mod = dest_time_mod.GetRasterBand(1).ReadAsArray()
            Array_time_myd = dest_time_myd.GetRasterBand(1).ReadAsArray()
            Array_lst_mod = dest_lst_mod.GetRasterBand(1).ReadAsArray()
            Array_lst_myd = dest_lst_myd.GetRasterBand(1).ReadAsArray()

            LST = Array_lst_mod
            Time = Array_time_mod
            VZA = Array_angle_mod

            LST = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_lst_myd, LST)
            Time = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_time_myd, Time)
            VZA = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_angle_myd, VZA)

            proj_ex = dest_angle_mod.GetProjection()
            geo_ex = dest_angle_mod.GetGeoTransform()

            PF.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)
            PF.Save_as_tiff(Time_file, Time, geo_ex, proj_ex)
            PF.Save_as_tiff(VZA_file, VZA, geo_ex, proj_ex)

        lst_files.append(LST_file)
        time_files.append(Time_file)
        vza_files.append(VZA_file)

    return (lst_files, time_files, vza_files)

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

def calc_dlat_dlon(geo_out, size_X, size_Y):
    """
    This functions calculated the distance between each pixel in meter.

    Parameters
    ----------
    geo_out: array
        geo transform function of the array
    size_X: int
        size of the X axis
    size_Y: int
        size of the Y axis

    Returns
    -------
    dlat: array
        Array containing the vertical distance between each pixel in meters
    dlon: array
        Array containing the horizontal distance between each pixel in meters
    """

    # Create the lat/lon rasters
    lon = np.arange(size_X + 1)*geo_out[1]+geo_out[0] - 0.5 * geo_out[1]
    lat = np.arange(size_Y + 1)*geo_out[5]+geo_out[3] - 0.5 * geo_out[5]

    dlat_2d = np.array([lat,]*int(np.size(lon,0))).transpose()
    dlon_2d =  np.array([lon,]*int(np.size(lat,0)))

    # Radius of the earth in meters
    R_earth = 6371000

    # Calculate the lat and lon in radians
    lonRad = dlon_2d * np.pi/180
    latRad = dlat_2d * np.pi/180

    # Calculate the difference in lat and lon
    lonRad_dif = abs(lonRad[:,1:] - lonRad[:,:-1])
    latRad_dif = abs(latRad[:-1] - latRad[1:])

    # Calculate the distance between the upper and lower pixel edge
    a = np.sin(latRad_dif[:,:-1]/2) * np.sin(latRad_dif[:,:-1]/2)
    clat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    dlat = R_earth * clat

    # Calculate the distance between the eastern and western pixel edge
    b = np.cos(latRad[1:,:-1]) * np.cos(latRad[:-1,:-1]) * np.sin(lonRad_dif[:-1,:]/2) * np.sin(lonRad_dif[:-1,:]/2)
    clon = 2 * np.arctan2(np.sqrt(b), np.sqrt(1-b))
    dlon = R_earth * clon

    return(dlat, dlon)