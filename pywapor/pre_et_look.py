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
import json
import xarray as xr
import glob
import warnings
import shutil

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1", 
        diagnostics = None, composite_length = "DEKAD"):

    print("\n#################")
    print("## PRE_ET_LOOK ##")
    print("#################")

#%%
    # Disable Warnings
    warnings.filterwarnings('ignore')
    all_files = dict()

    sdate = dat.strptime(startdate, "%Y-%m-%d").date()
    edate = dat.strptime(enddate, "%Y-%m-%d").date()
    epochs_info = create_dates(sdate, edate, composite_length)

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
    temp_folder = os.path.join(project_folder, "temporary")
    if diagnostics:
        diagnostics["folder"] = os.path.join(level_folder, "graphs")

    dl_args = (raw_folder, latlim, lonlim, startdate, enddate)

    unraw_file_templates = unraw_filepaths(epochs_info[1], level_folder, "{var}")

    #### NDVI ####
    print(f"\n#### NDVI ####")
    raw_ndvi_files = list()
    # Order is important! PROBV gets priority over MOD13, and MOD13 over MYD13.
    if "PROBAV" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.PROBAV.PROBAV_S5(*dl_args)[0])
        # raw_ndvi_files.append(glob.glob(os.path.join(project_folder, "RAW/PROBAV/NDVI/*.tif")))
    if "MOD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MOD13.NDVI(*dl_args))
    if "MYD13" in source_selection["NDVI"]:
        raw_ndvi_files.append(c.MYD13.NDVI(*dl_args))

    example_fh, example_ds, example_geoinfo = select_template(raw_ndvi_files)

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "ndvi",
        "var_unit": "-",
    }
    ds = g.compositer.main(cmeta, raw_ndvi_files, epochs_info, temp_folder, example_ds, 
                                diagnostics = diagnostics)

    all_files["ndvi"] = unraw_all(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

#%%

    #### ALBEDO ####
    print(f"\n#### ALBEDO ####")
    raw_albedo_files = list()
    # Order is important! PROBV gets priority over MDC43.
    if "PROBAV" in source_selection["ALBEDO"]:
        raw_albedo_files.append(c.PROBAV.PROBAV_S5(*dl_args)[1])
        # raw_albedo_files.append(glob.glob(r"/Volumes/Data/FAO/WaPOR_vs_pyWaPOR/pyWAPOR_long_test/RAW/PROBAV/Albedo/*.tif"))
    if "MDC43" in source_selection["ALBEDO"]:
        raw_albedo_files.append(c.MCD43.ALBEDO(*dl_args))

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "r0",
        "var_unit": "-",
    }
    ds = g.compositer.main(cmeta, raw_albedo_files, epochs_info, temp_folder, example_ds, diagnostics = diagnostics)
    all_files["r0"] = unraw_all(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

#%%
    #### PRECIPITATION ####
    print(f"\n#### PRECIPITATION ####")
    raw_precip_files = list()
    if "CHIRPS" in source_selection["PRECIPITATION"]:
        raw_precip_files.append(c.CHIRPS.daily(*dl_args))

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "nearest",
        "var_name": "P_24",
        "var_unit": "mm/day",
    }
    ds = g.compositer.main(cmeta, raw_precip_files, epochs_info, temp_folder, example_ds, diagnostics = diagnostics)
    all_files["P_24"] = unraw_all(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

#%%
    #### DEM ####
    print(f"\n#### DEM ####")
    if "SRTM" in source_selection["DEM"]:
        raw_dem_file = c.SRTM.DEM(*dl_args[:3])
    
    dem_file = unraw_filepaths(None, level_folder, "z", static = True)[0]

    cmeta = {
        "composite_type": None,
        "temporal_interp": None,
        "spatial_interp": "linear",
        "var_name": "z",
        "var_unit": "m",
    }
    ds = g.compositer.main(cmeta, [[raw_dem_file]], None, temp_folder, example_ds)
    all_files["z"] = unraw_all(cmeta["var_name"], ds, [dem_file], example_geoinfo)

    #### SLOPE ASPECT ####
    print(f"\n#### SLOPE ASPECT ####")
    slope_aspect(dem_file, level_folder, example_fh)

    #### LULC ####
    print(f"\n#### LULC ####")
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
            # print(year, raw_file, key)
            unraw_replace_values(raw_file, lulc_file_template.format(var = key, year = year), replace_values, example_fh)

#%%
    #### METEO ####
    print(f"\n#### METEO ####")
    if "MERRA2" in source_selection["METEO"]:
        meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
        raw_meteo_files = c.MERRA.daily_MERRA2(*dl_args, meteo_vars)
        raw_meteo_files = {**c.MERRA.daily_MERRA2(*dl_args, ['t2m'], data_type = ["min", "max"]), **raw_meteo_files}
        meteo_name_convertor = {'t2m-max': "t_air_max_24", 't2m-min': "t_air_min_24", 't2m': "t_air_24", 
                                'u2m': "u2m_24", 'v2m': "v2m_24",'q2m': "qv_24",'slp': "p_air_24_0"}
    elif "GEOS5" in source_selection["METEO"]:
        meteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']
        raw_meteo_files = c.GEOS.daily(*dl_args, meteo_vars)
        meteo_name_convertor = {'t2m': "t_air_24", 't2m-min': "t_air_min_24", 't2m-max': "t_air_max_24",  
                                'u2m': "u2m_24", 'v2m': "v2m_24",'slp': "p_air_24_0",'qv2m': "qv_24",}

    for var_name, raw_files in raw_meteo_files.items():

        print(f"## {var_name} ##")

        if "t2m" in var_name:
            raw_files = [g.lapse_rate.lapse_rate_temperature(x, dem_file) for x in raw_files if "_K_" in x]
        if var_name not in meteo_name_convertor.keys():
            continue
        units = {"t_air_max_24": "°C", "t_air_min_24":"°C", "t_air_24":"°C", 
          "u2m_24":"m/s", "v2m_24":"m/s","qv_24":"kg/kg","p_air_24_0":"kPa"}

        cmeta = {
            "composite_type": "mean",
            "temporal_interp": "linear",
            "spatial_interp": "linear",
            "var_name": meteo_name_convertor[var_name],
            "var_unit": units[meteo_name_convertor[var_name]],
        }
        ds = g.compositer.main(cmeta, [raw_files], epochs_info, temp_folder, example_ds, diagnostics = diagnostics)
        all_files[meteo_name_convertor[var_name]] = unraw_all(meteo_name_convertor[var_name], ds, unraw_file_templates, example_geoinfo)

#%%
    #### SE_ROOT ####
    ds_lst, ds_meteo, ds_ndvi, ds_temperature = pywapor.pre_se_root.main(project_folder, startdate, enddate, latlim, lonlim, level = source_selection)
    raw_se_root_files = pywapor.se_root.main(level_folder, ds_lst, ds_meteo, ds_ndvi, 
                                        ds_temperature, example_ds, example_geoinfo)

    cmeta = {
        "composite_type": 0.85,
        "temporal_interp": False,
        "spatial_interp": "nearest",
        "var_name": "se_root",
        "var_unit": "-",
    }
    ds = g.compositer.main(cmeta, [raw_se_root_files], epochs_info, temp_folder, example_ds, diagnostics = diagnostics)
    all_files["se_root"] = unraw_all(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

#%%
    #### LAT LON ####
    print(f"\n#### LAT LON ####")
    lat_lon(example_ds, level_folder, example_geoinfo)

    #### SOLAR RADIATION ####
    print(f"\n#### SOLAR RADIATION ####")
    if "MERRA2" in source_selection["TRANS"]:
        raw_ra24_files = c.MERRA.daily_MERRA2(*dl_args, ['swgnet'])['swgnet']

    cmeta = {
        "composite_type": "mean",
        "temporal_interp": "linear",
        "spatial_interp": "linear",
        "var_name": "ra_24",
        "var_unit": "-",
    }
    ds = g.compositer.main(cmeta, [raw_ra24_files], epochs_info, temp_folder, example_ds, diagnostics = diagnostics)
    all_files["ra_24"] = unraw_all(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)

    #### TEMP. AMPLITUDE #### # TODO add source selection and remove year data
    print(f"\n#### TEMP. AMPLITUDE ####")
    raw_temp_ampl_file = os.path.join(raw_folder, "GLDAS", "Temp_Amplitudes_global.tif")
    download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", raw_temp_ampl_file)
    temp_ampl_file_template = unraw_filepaths(None, level_folder, "t_amp_year_{year}", static = True)[0]
    raw_temp_ampl_files = [(year, raw_temp_ampl_file) for year in range(sdate.year, edate.year + 1)]
    for year, raw_file in raw_temp_ampl_files:
        unraw(raw_file, temp_ampl_file_template.format(year = year), example_fh, 6)

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

    for fh in glob.glob(os.path.join(temp_folder, "*.nc")):
        os.remove(fh)
    if len(os.listdir(temp_folder)) == 0:
        shutil.rmtree(temp_folder)

    os.chdir(project_folder)

    print("\n#################")
    print("## PRE_ET_LOOK ##")
    print("##### DONE ######\n")
#%%
    return all_files

def unraw_all(var_name, ds, unraw_file_templates, example_geoinfo):
    unrawed_files = list()
    assert ds.epoch.size == len(unraw_file_templates)
    for i, fh in enumerate(unraw_file_templates):
        if not -9999 in ds["epoch"].values:
            fh_date = pd.Timestamp(dat.strptime(fh.split("_")[-1], "%Y%m%d.tif"))
            assert fh_date == pd.Timestamp(ds.epoch_starts.isel(epoch = i).values)
        array = np.copy(ds.composite.isel(epoch = i).values)
        real_fh = fh.format(var = var_name)
        PF.Save_as_tiff(real_fh, array, example_geoinfo[0], example_geoinfo[1])
        unrawed_files.append(real_fh)
    ds.close()
    ds = None
    return unrawed_files

def select_template(fhs):
    fhs = [val for sublist in fhs for val in sublist]

    sizes = [gdal.Open(fh).RasterXSize * gdal.Open(fh).RasterYSize for fh in fhs]
    idx = np.argmax(sizes)

    example_fh = fhs[idx]
    example_ds = xr.open_dataset(example_fh).isel(band = 0).drop_vars(["band", "spatial_ref"]).rename({"x": "lon", "y": "lat"})
    example_geoinfo = PF.get_geoinfo(example_fh)

    return example_fh, example_ds, example_geoinfo

def slope_aspect(dem_file, project_folder, template_file):

    slope_file = unraw_filepaths(None, project_folder, "slope_deg", static = True)[0]
    aspect_file = unraw_filepaths(None, project_folder, "aspect_deg", static = True)[0]

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

def lat_lon(ds, project_folder, example_geoinfo):

    lat_file = unraw_filepaths(None, project_folder, "lat_deg", static = True)[0]
    lon_file = unraw_filepaths(None, project_folder, "lon_deg", static = True)[0]

    lat_deg = np.rot90(np.tile(ds.lat, ds.lon.size).reshape(ds.band_data.shape[::-1]), 3)
    lon_deg = np.tile(ds.lon, ds.lat.size).reshape(ds.band_data.shape)
    if not os.path.exists(lon_file):
        PF.Save_as_tiff(lon_file, lon_deg, example_geoinfo[0], example_geoinfo[1])
    if not os.path.exists(lat_file):
        PF.Save_as_tiff(lat_file, lat_deg, example_geoinfo[0], example_geoinfo[1])

    return lat_file, lon_file

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
    periods_end = dates[1:][periods_start < np.datetime64(edate)]
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

def download_file_from_google_drive(id, destination):
    
    if os.path.isfile(destination):
        return

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

    return           

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

# if __name__ == "__main__":

    # project_folder = r"/Volumes/Data/FAO/WaPOR_vs_pyWaPOR/pyWAPOR_v1"
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    # startdate = "2021-07-01"
    # enddate = "2021-07-10"
    # composite_length = "DEKAD"

    # level = {
    #     "METEO": ["GEOS5"],
    #     "NDVI": ["MOD13", "MYD13", "PROBAV"],
    #     "ALBEDO": ["PROBAV"],
    #     "LST": ["MOD11", "MYD11"],
    #     "LULC": ["WAPOR"],
    #     "DEM": ["SRTM"],
    #     "PRECIPITATION": ["CHIRPS"],
    #     "TRANS": ["MERRA2"],
    # }
    # level = "level_2"

    # diagnostics = { # label          # lat      # lon
    #                 "water":	    (29.44977,	30.58215),
    #                 "desert":	    (29.12343,	30.51222),
    #                 "agriculture":	(29.32301,	30.77599),
    #                 "urban":	    (29.30962,	30.84109),
    #                 }
    # diagnostics = None

    # all_files = main(project_folder, startdate, enddate, latlim, lonlim, level = level, 
    #     diagnostics = diagnostics, composite_length = composite_length)