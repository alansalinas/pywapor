# -*- coding: utf-8 -*-
"""
"""
import os
from datetime import datetime as dat
from osgeo import gdal
import pandas as pd
import numpy as np
import requests
import pywapor
import pywapor.general.processing_functions as PF
import pywapor.collect as c
import pywapor.general as g
from pathlib import Path
import json

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1"):

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
    level_folder = os.path.join(project_folder, level)

    dl_args = (raw_folder, latlim, lonlim, startdate, enddate)

    unraw_file_templates = unraw_filepaths(startdate, enddate, level_folder, "{var}")

    #### NDVI #### 
    raw_ndvi_files = list()
    # Order is important! PROBV gets priority over MOD13, and MOD13 over MYD13.
    if "PROBAV" in source_selection["NDVI"]:
        raw_ndvi_files += c.PROBAV.PROBAV_S5(*dl_args)[0]
    if "MOD13" in source_selection["NDVI"]:
        raw_ndvi_files += c.MOD13.NDVI(*dl_args)
    if "MYD13" in source_selection["NDVI"]:
        raw_ndvi_files += c.MYD13.NDVI(*dl_args)

    template_file = select_template(raw_ndvi_files)

    unraw_all("NDVI", unraw_file_templates, raw_ndvi_files, template_file, 1)

    #### ALBEDO ####
    raw_albedo_files = list()
    # Order is important! PROBV gets priority over MDC43.
    if "PROBAV" in source_selection["ALBEDO"]:
        raw_albedo_files += c.PROBAV.PROBAV_S5(*dl_args)[1]
    if "MDC43" in source_selection["ALBEDO"]:
        raw_albedo_files += c.MCD43.ALBEDO(*dl_args)

    unraw_all("ALBEDO", unraw_file_templates, raw_albedo_files, template_file, 1)

    #### LST ####
    raw_lst_files = list()
    if "MOD11" in source_selection["LST"]:
        raw_lst_files.append(c.MOD11.LST(*dl_args))
    if "MYD11" in source_selection["LST"]:
        raw_lst_files.append(c.MYD11.LST(*dl_args))
    
    raw_lst_files, raw_time_files = combine_lst(raw_lst_files)

    unraw_all("LST", unraw_file_templates, raw_lst_files, template_file, 2)

    #### TIME ####
    time_files = unraw_all("Time", unraw_file_templates, raw_time_files, template_file, 1)

    #### PRECIPITATION ####
    if "CHIRPS" in source_selection["PRECIPITATION"]:
        raw_precip_files = c.CHIRPS.daily(*dl_args)

    unraw_all("Precipitation", unraw_file_templates, raw_precip_files, template_file, 6)

    #### DEM ####
    if "SRTM" in source_selection["DEM"]:
        raw_dem_file = c.SRTM.DEM(*dl_args[:3])
    
    dem_file = unraw_filepaths("", "", level_folder, "DEM", static = True)[0]
    unraw(raw_dem_file, dem_file, template_file, 4)

    #### SLOPE ASPECT ####
    slope_aspect(dem_file, level_folder, template_file)

    #### LULC ####
    if "GLOBCOVER" in source_selection["LULC"]:
        raw_lulc_file = c.Globcover.Landuse(*dl_args[:3])
        raw_lulc_files = [(year, raw_lulc_file) for year in range(sdate.year, edate.year + 1)]
    elif "WAPOR" in source_selection["LULC"]:
        raw_lulc_files = c.WAPOR.Get_Layer(*dl_args[:3], sdate.strftime("%Y-01-01"), edate.strftime("%Y-12-31"), 'L1_LCC_A')
        raw_lulc_files = [(dat.strptime(os.path.split(fp)[-1], "L1_LCC_A_WAPOR_YEAR_%Y.%m.%d.tif").year, fp) for fp in raw_lulc_files]

    lulc_values = g.landcover_converter.get_lulc_values()
    
    lulc_file_template = unraw_filepaths("", "", level_folder, "{var}_{year}", static = True)[0]
    for year, raw_file in raw_lulc_files:
        for key, replace_values in lulc_values[source_selection["LULC"][0]].items():
            print(year, raw_file, key)
            unraw_replace_values(raw_file, lulc_file_template.format(var = key, year = year), replace_values, template_file)

    #### METEO ####
    if "MERRA2" in source_selection["METEO"]:
        freq = "H"
        periods = calc_periods(time_files, freq)
        meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
        c.MERRA.daily_MERRA2(*dl_args, meteo_vars)
        c.MERRA.daily_MERRA2(*dl_args, ['t2m'], data_type = ["mean", "min", "max"])
        for sd, ed, period in periods:
            c.MERRA.hourly_MERRA2(*dl_args[:3], sd, ed, meteo_vars, [int(period)])
    elif "GEOS5" in source_selection["METEO"]:
        freq = "3H"
        periods = calc_periods(time_files, freq)
        meteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']
        c.GEOS.daily(*dl_args, meteo_vars)
        for sd, ed, period in periods:
            c.GEOS.three_hourly(*dl_args[:3], sd, ed, meteo_vars, [int(period)])

    meteo_files_template = unraw_filepaths(startdate, enddate, level_folder, "{var}")
    
    raw_meteo_paths = g.variables.get_raw_meteo_paths()

    for meteo_file_template in meteo_files_template:
        date = dat.strptime(os.path.split(meteo_file_template)[-1], "{var}_%Y%m%d.tif")
        date_str = date.strftime("%Y.%m.%d")
        period = [x[2] for x in periods if x[0] == date][0] 
        hour = int((period - 1) * {"3H": 3, "H": 1}[freq])
        hour_str = str(hour).zfill(2)
        for key, path in raw_meteo_paths[source_selection["METEO"][0]].items():
            unraw_file = meteo_file_template.format(var = key)
            if "wind" in key:
                raw_file_u = os.path.join(*path[0]).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                raw_file_v = os.path.join(*path[1]).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                u_wind = PF.reproj_file(raw_file_u, template_file, 6)
                v_wind = PF.reproj_file(raw_file_v, template_file, 6)
                wind = np.sqrt(u_wind**2 + v_wind**2)
                geo_ex, proj_ex = PF.get_geoinfo(template_file)[0:2]
                PF.Save_as_tiff(unraw_file, wind, geo_ex, proj_ex)
            elif "tair" in key:
                raw_file = os.path.join(*path).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                tair = lapse_rate_temp(raw_file, dem_file)
                geo_ex, proj_ex = PF.get_geoinfo(template_file)[0:2]
                PF.Save_as_tiff(unraw_file, tair, geo_ex, proj_ex)
            else:
                raw_file = os.path.join(*path).format(raw_folder = raw_folder, date = date_str, hour = hour_str)
                unraw(raw_file, unraw_file, template_file, method = 6)

    #### LAT LON ####
    lat_file = lat_lon(level_folder, template_file)[0]

    #### TRANS ####
    if "MERRA2" in source_selection["TRANS"]:
        raw_trans_files = c.MERRA.daily_MERRA2(*dl_args, ['swgnet'])[0]
    
    trans_files = unraw_filepaths(startdate, enddate, level_folder, "Trans_24")
    trans_files = [(int(dat.strptime(os.path.split(fp)[-1], "Trans_24_%Y%m%d.tif").strftime("%j")), fp) for fp in trans_files]

    for (doy, unraw_file), raw_file in zip(trans_files, raw_trans_files):
        ra24_flat = calc_ra24_flat(lat_file, doy)
        array = PF.reproj_file(raw_file, template_file, 6) / ra24_flat
        geo_ex, proj_ex = PF.get_geoinfo(template_file)[0:2]
        PF.Save_as_tiff(unraw_file, array, geo_ex, proj_ex)

    #### TEMP. AMPLITUDE ####
    raw_temp_ampl_file = os.path.join(raw_folder, "GLDAS", "Temp_Amplitudes_global.tif")
    download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", raw_temp_ampl_file)
    temp_ampl_file_template = unraw_filepaths("", "", level_folder, "Tair_amp_{year}", static = True)[0]
    raw_temp_ampl_files = [(year, raw_temp_ampl_file) for year in range(sdate.year, edate.year + 1)]
    for year, raw_file in raw_temp_ampl_files:
        unraw(raw_file, temp_ampl_file_template.format(year = year), template_file, 6)

    #### METADATA ####
    metadata = dict()
    metadata["pywapor_version"] = pywapor.__version__
    metadata["created"] = dat.now().strftime("%m/%d/%Y, %H:%M:%S")
    metadata["template_file"] = template_file
    metadata["geotransform"] = "[{0}, {1}, {2}, {3}, {4}, {5}]".format(*PF.get_geoinfo(template_file)[0])
    metadata["resolution"] = "[{0}, {1}]".format(*PF.get_geoinfo(template_file)[2:])
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

def match_to_nearest(raw_files, files):
    raw_dates = [dat.strptime(os.path.split(fp)[-1].split("_")[-1], "%Y.%m.%d.tif") for fp in raw_files]
    dates = [dat.strptime(os.path.split(fp)[-1].split("_")[-1], "%Y%m%d.tif") for fp in files]
    find_idx = lambda d: np.argmin([np.abs((raw_date - d).days) for raw_date in raw_dates])
    matched_raw_ndvi_files = [raw_files[find_idx(d)] for d in dates]
    return matched_raw_ndvi_files

def select_template(fhs):
    sizes = [gdal.Open(fh).RasterXSize * gdal.Open(fh).RasterYSize for fh in fhs]
    idx = np.argmax(sizes)
    return fhs[idx]
    
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

    slope_file = unraw_filepaths("", "", project_folder, "Slope", static = True)[0]
    aspect_file = unraw_filepaths("", "", project_folder, "Aspect", static = True)[0]

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

    lat_file = unraw_filepaths("", "", project_folder, "Lat", static = True)[0]
    lon_file = unraw_filepaths("", "", project_folder, "Lon", static = True)[0]

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

def unraw_filepaths(startdate, enddate, project_folder, var, static = False):
    filepaths = list()
    if not static:
        dates = pd.date_range(startdate, enddate, freq="D")
        for date in dates:
            date_str = date.strftime("%Y%m%d")
            fp = os.path.join(project_folder,
                            date_str, f"{var}_{date_str}.tif")
            filepaths.append(fp)
    else:
        fp = os.path.join(project_folder,
                        "static", f"{var}.tif")
        filepaths.append(fp)   
    return filepaths

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

    tdown = pywapor.et_look_v2.meteo.disaggregate_air_temperature(tempe, dem_down, dem_up_ave)
    # _plot_array(tdown)

    # test_tair = r"/Volumes/Data/pre_et_look_ORIGINAL/ETLook_input_MODIS/20190706/tair_24_20190706.tif"
    # tempe_test = PF.open_as_array(test_tair)
    # _plot_array(tempe_test - tdown)

    return tdown

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
