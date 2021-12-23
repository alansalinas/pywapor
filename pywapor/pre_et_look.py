# -*- coding: utf-8 -*-
"""
"""
import os
import json
import glob
import warnings
import shutil
import requests
import pywapor
import xarray as xr
import pandas as pd
import numpy as np
import pywapor.general as g
import pywapor.enhancers as enhance
import pywapor.general.processing_functions as pf
import pywapor.general.pre_defaults as defaults
from datetime import datetime as dat
from osgeo import gdal
from pywapor.collect.downloader import download_sources
from pywapor.general.logger import log, adjust_logger

def main(project_folder, startdate, enddate, latlim, lonlim, level = "level_1", 
        diagnostics = None, composite_length = "DEKAD", extra_sources = None,
        extra_source_locations = None):

    if not os.path.exists(project_folder):
        os.makedirs(project_folder)

    log_write = True
    log_level = "INFO"
    adjust_logger(log_write, project_folder, log_level)

    log.info("> PRE_ET_LOOK").add()

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
    succes = g.tests.check_source_selection(source_selection, sdate, edate, extra_sources)[1]
    assert succes, "invalid source_selection"

    raw_folder = os.path.join(project_folder, "RAW")
    level_folder = os.path.join(project_folder, level)
    temp_folder = os.path.join(project_folder, "temporary")
    if isinstance(diagnostics, dict):
        diagnostics["folder"] = os.path.join(level_folder, "graphs")

    dl_args = {
                "Dir": raw_folder, 
                "latlim": latlim, 
                "lonlim": lonlim, 
                "Startdate": startdate, 
                "Enddate": enddate,
                }

    unraw_file_templates = unraw_filepaths(epochs_info[1], level_folder, "{var}")

    #### NDVI ####
    var = "NDVI"
    log.info(f"# {var}")
    raw_ndvi_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    example_fh, example_ds, example_geoinfo = select_template(raw_ndvi_files)
    all_files[var] = unraw_all(var, raw_ndvi_files, epochs_info, temp_folder, 
                example_ds, unraw_file_templates, example_geoinfo)

    #### ALBEDO ####
    var = "ALBEDO"
    log.info(f"# {var}")
    raw_albedo_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    all_files[var] = unraw_all(var, raw_albedo_files, epochs_info, temp_folder, 
                example_ds, unraw_file_templates, example_geoinfo)

    #### PRECIPITATION ####
    var = "PRECIPITATION"
    log.info(f"# {var}")
    raw_precip_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    all_files[var] = unraw_all(var, raw_precip_files, epochs_info, temp_folder, 
                example_ds, unraw_file_templates, example_geoinfo)

    #### DEM ####
    var = "DEM"
    log.info(f"# {var}")
    raw_dem_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    dem_file = unraw_filepaths(None, level_folder, "z", static = True)
    all_files[var] = unraw_all(var, raw_dem_files, epochs_info, temp_folder, 
                example_ds, dem_file, example_geoinfo)

    #### LULC ####
    var = "LULC"
    log.info(f"# {var}")
    raw_lulc_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    # lulc_values = enhance.landcover_converter.get_lulc_values()
    all_files[var] = unraw_all(var, raw_lulc_files, epochs_info, temp_folder, 
                example_ds, unraw_file_templates, example_geoinfo)

    #### METEO ####
    log.info("> METEO").add()
    meteo_vars = ['t_air_24', 't_air_min_24', 't_air_max_24', 
                  'u2m_24', 'v2m_24', 'p_air_24_0', 'qv_24']
    for meteo_var in meteo_vars:
        raw_meteo_files = download_sources(meteo_var, source_selection["METEO"], 
                                            dl_args, extra_source_locations)
        all_files[meteo_var] = unraw_all(meteo_var, raw_meteo_files, epochs_info, 
                                    temp_folder, example_ds, unraw_file_templates, 
                                    example_geoinfo)
    log.sub().info("< METEO")

    #### SE_ROOT ####
    ds_lst, ds_meteo, ds_ndvi, ds_temperature = pywapor.pre_se_root.main(project_folder, startdate, enddate, latlim, 
                                                                        lonlim, level = source_selection, extra_sources = extra_sources,
                                                                        extra_source_locations = extra_source_locations)
    raw_se_root_files = pywapor.se_root.main(level_folder, ds_lst, ds_meteo, ds_ndvi, 
                                        ds_temperature, example_ds, example_geoinfo)

    var = "SE_ROOT"
    all_files[var] = unraw_all(var, raw_se_root_files, epochs_info, temp_folder, 
                                example_ds, unraw_file_templates, example_geoinfo)

    #### SOLAR RADIATION ####
    var = "SOLAR_RADIATION"
    log.info(f"# {var}")
    raw_ra24_files = download_sources(var, source_selection[var], 
                                        dl_args, extra_source_locations)
    all_files[var] = unraw_all(var, raw_ra24_files, epochs_info, temp_folder, 
                                example_ds, unraw_file_templates, example_geoinfo)

    #### SLOPE ASPECT ####
    log.info("# SLOPE ASPECT")
    slope_aspect(dem_file, level_folder, example_fh)

    #### LAT LON ####
    log.info("# LAT LON")
    lat_lon(example_ds, level_folder, example_geoinfo)

    #### TEMP. AMPLITUDE #### # TODO add source selection and remove year data
    log.info("# TEMP.-AMPLITUDE")
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
    metadata["geotransform"] = "[{0}, {1}, {2}, {3}, {4}, {5}]".format(*pf.get_geoinfo(example_fh)[0])
    metadata["resolution"] = "[{0}, {1}]".format(*pf.get_geoinfo(example_fh)[2:])
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

    log.sub().info("< PRE_ET_LOOK")

    return all_files

def unraw_all(var, raw_files, epochs_info, temp_folder, example_ds, unraw_file_templates, example_geoinfo):
    cmeta = defaults.composite_defaults()[var]
    ds = g.compositer.main(cmeta, raw_files, epochs_info, temp_folder, example_ds)
    if isinstance(diagnostics, dict):
        ds_diags = g.compositer.main(cmeta, raw_files, epochs_info, None, example_ds, 
                                    lean_output = False, diagnostics = diagnostics)
    files = ds_to_geotiff(cmeta["var_name"], ds, unraw_file_templates, example_geoinfo)
    return files

def ds_to_geotiff(var_name, ds, unraw_file_templates, example_geoinfo):
    unrawed_files = list()
    assert ds.epoch.size == len(unraw_file_templates)
    for i, fh in enumerate(unraw_file_templates):
        if not -9999 in ds["epoch"].values:
            fh_date = pd.Timestamp(dat.strptime(fh.split("_")[-1], "%Y%m%d.tif"))
            assert fh_date == pd.Timestamp(ds.epoch_starts.isel(epoch = i).values)
        array = np.copy(ds.composite.isel(epoch = i).values)
        real_fh = fh.format(var = var_name)
        pf.Save_as_tiff(real_fh, array, example_geoinfo[0], example_geoinfo[1])
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
    example_geoinfo = pf.get_geoinfo(example_fh)

    resolution = pywapor.et_look.get_geoinfo(example_fh)[0]
    log.info(f"--> Resampling resolution is ~{resolution:.0f} meter.")

    return example_fh, example_ds, example_geoinfo

def slope_aspect(dem_file, project_folder, template_file):

    slope_file = unraw_filepaths(None, project_folder, "slope_deg", static = True)[0]
    aspect_file = unraw_filepaths(None, project_folder, "aspect_deg", static = True)[0]

    if not os.path.exists(slope_file) or not os.path.exists(aspect_file):

        dem = pf.open_as_array(dem_file)

        # constants
        geo_ex, proj_ex, size_x_ex, size_y_ex = pf.get_geoinfo(template_file)
        dlat, dlon = pf.calc_dlat_dlon(geo_ex, size_x_ex, size_y_ex)            

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
        pf.Save_as_tiff(slope_file, slope, geo_ex, proj_ex)
        pf.Save_as_tiff(aspect_file, aspect, geo_ex, proj_ex)

    return slope_file, aspect_file

def lat_lon(ds, project_folder, example_geoinfo):

    lat_file = unraw_filepaths(None, project_folder, "lat_deg", static = True)[0]
    lon_file = unraw_filepaths(None, project_folder, "lon_deg", static = True)[0]

    lat_deg = np.rot90(np.tile(ds.lat, ds.lon.size).reshape(ds.band_data.shape[::-1]), 3)
    lon_deg = np.tile(ds.lon, ds.lat.size).reshape(ds.band_data.shape)
    if not os.path.exists(lon_file):
        pf.Save_as_tiff(lon_file, lon_deg, example_geoinfo[0], example_geoinfo[1])
    if not os.path.exists(lat_file):
        pf.Save_as_tiff(lat_file, lat_deg, example_geoinfo[0], example_geoinfo[1])

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

def unraw(raw_file, unraw_file, template_file, method):
    if not os.path.exists(unraw_file) and os.path.exists(raw_file):
        geo_ex, proj_ex = pf.get_geoinfo(template_file)[0:2]
        array = pf.reproj_file(raw_file, template_file, method)
        pf.Save_as_tiff(unraw_file, array, geo_ex, proj_ex)

def unraw_replace_values(raw_file, unraw_file, replace_values, template_file):
    if not os.path.exists(unraw_file) and os.path.exists(raw_file):
        array = pf.reproj_file(raw_file, template_file, 1)
        replaced_array = np.ones_like(array) * np.nan
        for key, value in replace_values.items():
            replaced_array[array == key] = value
        geo_ex, proj_ex = pf.get_geoinfo(template_file)[0:2]
        pf.Save_as_tiff(unraw_file, replaced_array, geo_ex, proj_ex)

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

def check_extra_product_names(extra_products):
    """Product names need to be unique and can't contain underscores.
    Variables for which new products can be added are 'NDVI', 'ALBEDO' and 
    'LST'.

    Parameters
    ----------
    extra_products : dict
        Keys are variables (e.g. 'NDVI'), values are lists of new products.

    """
    valid_keys = ["NDVI", "ALBEDO", "LST"]

    assert np.all([x in valid_keys for x in extra_products.keys()])

    all_products = np.array(list(extra_products.values())).flatten()

    assert all_products.size == np.unique(all_products).size
    assert np.all(["_" not in x for x in all_products])

if __name__ == "__main__":

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"
    composite_length = "DEKAD"
    level = "level_1"
    extra_sources = None
    extra_source_locations = None

    diagnostics = { # label          # lat      # lon
                    "water":	    (29.44977,	30.58215),
                    "desert":	    (29.12343,	30.51222),
                    "agriculture":	(29.32301,	30.77599),
                    "urban":	    (29.30962,	30.84109),
                    }
    # diagnostics = None

    level = {'METEO': ['MERRA2'],
            'NDVI': ['MYD13'],
            'ALBEDO': ['MCD43'],
            'LST': ['MOD11', 'MYD11'],
            'LULC': ['WAPOR'],
            'DEM': ['SRTM'],
            'PRECIPITATION': ['CHIRPS'],
            'SOLAR_RADIATION': ['MERRA2'],
            "name": "level_1"}

    # # Define new products.
    # extra_sources =  {"NDVI":      ["LS7NDVI", "LS8NDVI"],
    #                     "LST":     ["LS7LST", "LS8LST"],
    #                     "ALBEDO":  ["LS7ALBEDO", "LS8ALBEDO"]}

    # # Give product folder.
    # extra_source_locations = {
    #     "LS7NDVI": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/NDVI",
    #     "LS8NDVI": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/NDVI",
    #     "LS7LST": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LST",
    #     "LS8LST": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LST",
    #     "LS7ALBEDO": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/ALBEDO",
    #     "LS8ALBEDO": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/ALBEDO",
    # }

    all_files = main(project_folder, startdate, enddate, latlim, lonlim, level = level, 
        diagnostics = diagnostics, composite_length = composite_length, extra_sources = extra_sources,
        extra_source_locations = extra_source_locations)

    # all_final_files = pywapor.et_look.main(project_folder, startdate)

    # param = "SOLAR_RADIATION"
    # sources = ["GEOS5"]

    # dl_args = {
    #             "Dir": os.path.join(project_folder, "RAW"), 
    #             "latlim": latlim, 
    #             "lonlim": lonlim, 
    #             "Startdate": startdate, 
    #             "Enddate": enddate,
    #             }

    # files = download_sources(param, sources, dl_args)



