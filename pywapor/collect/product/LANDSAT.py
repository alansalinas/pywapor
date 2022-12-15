"""
LXSS_LLLL_PPPRRR_YYYYMMDD_yyyymmdd_CX_TX
(e.g., LC08_L2SP_039037_20150728_20200318_02_T1)
L Landsat
X Sensor (“O” = OLI; “T” = TIRS; “C” = OLI/TIRS)
SS Satellite (“08” = Landsat 8, “09” = Landsat 9)
LLLL Processing correction level (“L2SP” if SR and ST are generated or “L2SR”
if ST could not be generated) PPP Path
RRR Row
YYYY Year of acquisition
MM Month of acquisition
DD Day of acquisition
yyyy Year of Level 2 processing
mm Month of Level 2 processing
dd Day of Level 2 processing
CX Collection number (“02”)
TX Collection category ( “T1” = Tier 1; “T2” = Tier 2)
"""

import requests as r
import pywapor
import tarfile
import re
import shutil
import glob
import json
import warnings
import os
from functools import partial
import numpy as np
import xarray as xr
import datetime
import rasterio
from pywapor.general.logger import log
from pywapor.general.processing_functions import open_ds, save_ds, remove_ds, adjust_timelim_dtype, make_example_ds
from pywapor.collect.protocol.crawler import download_urls
from pywapor.enhancers.apply_enhancers import apply_enhancer
import xmltodict
from pywapor.general import bitmasks

def calc_normalized_difference(ds, var, bands = ["nir", "red"]):
    """Calculate the normalized difference of two bands.

    Parameters
    ----------
    ds : xr.Dataset
        Input data.
    var : str
        Name of the variable in which to store the normalized difference.
    bands : list, optional
        The two bands to use to calculate the norm. difference, by default ["nir", "red"].

    Returns
    -------
    xr.Dataset
        Output data.
    """
    if np.all([x in ds.data_vars for x in bands]):
        da = (ds[bands[0]] - ds[bands[1]]) / (ds[bands[0]] + ds[bands[1]])
        ds[var] = da.clip(-1, 1)
    else:
        log.warning(f"--> Couldn't calculate `{var}`, `{'` and `'.join([x for x in bands if x not in ds.data_vars])}` is missing.")
    return ds

def mask_uncertainty(ds, var, max_uncertainty = 3):
    if "lst_qa" in ds.data_vars:
        ds[var] = ds[var].where(ds["lst_qa"] <= max_uncertainty, np.nan)
    else:
        log.warning(f"--> Couldn't mask uncertain values.")
    return ds

def mask_invalid(ds, var):
    if "valid_range" in ds[var].attrs.keys():
        ds[var] = xr.where(
                            (ds[var] >= ds[var].valid_range[0]) & 
                            (ds[var] <= ds[var].valid_range[1]) &
                            (ds[var] != ds[var]._FillValue), 
                            ds[var], np.nan, keep_attrs=True)
    else:
        log.warning(f"--> Couldn't mask invalid values, since `valid_range` is not defined.")
    return ds

def apply_qa(ds, var, pixel_qa_flags = None, radsat_qa_flags = None):
    
    masks = list()

    if ("pixel_qa" in ds.data_vars) and not isinstance(pixel_qa_flags, type(None)):
        pixel_qa_bits = bitmasks.get_pixel_qa_bits(2, int(product_name[-1]), 2)
        mask1 = bitmasks.get_mask(ds["pixel_qa"], pixel_qa_flags, pixel_qa_bits)
        masks.append(mask1)

    if ("radsat_qa" in ds.data_vars) and not isinstance(radsat_qa_flags, type(None)):
        radsat_qa_bits = bitmasks.get_radsat_qa_bits(2, int(product_name[-1]), 2)
        radsat_qa_flags = list(radsat_qa_bits.keys())
        mask2 = bitmasks.get_mask(ds["radsat_qa"], radsat_qa_flags, radsat_qa_bits)
        masks.append(mask2)

    if len(masks) >= 1:
        mask = np.invert(np.any(masks, axis = 0))
        ds[var] = ds[var].where(mask)

    return ds

def scale_data(ds, var):
    scale = getattr(ds[var], "scale_factor", 1.0)
    offset = getattr(ds[var], "add_offset", 0.0)
    ds[var] = ds[var] * scale + offset
    return ds

def default_vars(product_name, req_vars):
    # {x: [ds[x].dims, ds[x].long_name] for x in ds.data_vars}

    pixel_qa_flags_89 = ["dilated_cloud", "cirrus", "cloud", "cloud_shadow", "snow"]
    pixel_qa_flags_457 = ["dilated_cloud", "cloud", "cloud_shadow", "snow"]

    variables = {

        "LC08": {
            'sr_band1':     [('YDim_sr_band1', 'XDim_sr_band1'), 'coastal',[mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band1"])]],
            'sr_band2':     [('YDim_sr_band2', 'XDim_sr_band2'), 'blue', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band2"])]],
            'sr_band3':     [('YDim_sr_band3', 'XDim_sr_band3'), 'green', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band3"])]],
            'sr_band4':     [('YDim_sr_band4', 'XDim_sr_band4'), 'red', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band4"])]],
            'sr_band5':     [('YDim_sr_band5', 'XDim_sr_band5'), 'nir', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band5"])]],
            'sr_band6':     [('YDim_sr_band6', 'XDim_sr_band6'), 'swir1', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band6"])]],
            'sr_band7':     [('YDim_sr_band7', 'XDim_sr_band7'), 'swir2', [mask_invalid, scale_data, partial(apply_qa, pixel_qa_flags = pixel_qa_flags_89, radsat_qa_flags = ["terrain_occlusion", "saturated_band7"])]],
            'st_band10':    [('YDim_st_band10', 'XDim_st_band10'), 'lst', [mask_invalid, scale_data]],

            # 'st_trad':      [('YDim_st_trad', 'XDim_st_trad'), 'thermal_radiance'],
            # 'st_urad':      [('YDim_st_urad', 'XDim_st_urad'), 'upwell_radiance'],
            # 'st_drad':      [('YDim_st_drad', 'XDim_st_drad'), 'downwell_radiance'],
            # 'st_atran':     [('YDim_st_atran', 'XDim_st_atran'), 'atmospheric_transmittance'],
            # 'st_emis':      [('YDim_st_emis', 'XDim_st_emis'), 'emissivity_estimated_from_ASTER_GED'],
            # 'st_emsd':      [('YDim_st_emsd', 'XDim_st_emsd'), 'emissivity_standard_deviation'],
            # 'st_cdist':     [('YDim_st_cdist', 'XDim_st_cdist'), 'cloud_distance'],
            # 'sr_aerosol':   [('YDim_sr_aerosol', 'XDim_sr_aerosol'), 'surface_reflectance_aerosol_mask'],

            'st_qa':        [('YDim_st_qa', 'XDim_st_qa'), 'lst_qa', [mask_invalid, scale_data]],
            'qa_pixel':     [('YDim_qa_pixel', 'XDim_qa_pixel'), 'pixel_qa', []],
            'qa_radsat':    [('YDim_qa_radsat', 'XDim_qa_radsat'), 'radsat_qa', []],
        },

    }

    req_dl_vars = {

        "LC08": {
            'coastal': ['sr_band1', 'qa_pixel', 'qa_radsat'],
            'blue': ['sr_band2', 'qa_pixel', 'qa_radsat'],
            'green': ['sr_band3', 'qa_pixel', 'qa_radsat'],
            'red': ['sr_band4', 'qa_pixel', 'qa_radsat'],
            'nir': ['sr_band5', 'qa_pixel', 'qa_radsat'],
            'swir1': ['sr_band6', 'qa_pixel', 'qa_radsat'],
            'swir2': ['sr_band7', 'qa_pixel', 'qa_radsat'],
            'lst': ['st_band10', 'st_qa'],

            # 'thermal_radiance': ['st_trad'],
            # 'upwell_radiance': ['st_urad'],
            # 'downwell_radiance': ['st_drad'],
            # 'atmospheric_transmittance': ['st_atran'],
            # 'emissivity_estimated_from_ASTER_GED': ['st_emis'],
            # 'emissivity_standard_deviation': ['st_emsd'],
            # 'cloud_distance': ['st_cdist'],
            # 'surface_reflectance_aerosol_mask': ['sr_aerosol'],

            'lst_qa': ['st_qa'],
            'pixel_qa': ['qa_pixel'],
            'radsat_qa': ['qa_radsat'],

            'ndvi': ['sr_band4', 'sr_band5', 'qa_pixel', 'qa_radsat']
        },

    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars):

    post_processors = {
        "LC08": {
            'coastal': [],
            'blue': [],
            'green': [],
            'red': [],
            'nir': [],
            'swir1': [],
            'swir2': [],
            'lst': [mask_uncertainty],

            # 'thermal_radiance': [],
            # 'upwell_radiance': [],
            # 'downwell_radiance': [],
            # 'atmospheric_transmittance': [],
            # 'emissivity_estimated_from_ASTER_GED': [],
            # 'emissivity_standard_deviation': [],
            # 'cloud_distance': [],
            # 'surface_reflectance_aerosol_mask': [],

            'lst_qa': [],
            'pixel_qa': [],
            'radsat_qa': [],

            'ndvi': [calc_normalized_difference]
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def espa_api(endpoint, verb='get', body=None, uauth=None):
    """ Suggested simple way to interact with the ESPA JSON REST API """
    # auth_tup = uauth if uauth else (username, password)
    host = 'https://espa.cr.usgs.gov/api/v1/'
    response = getattr(r, verb)(host + endpoint, auth=uauth, json=body)
    print('{} {}'.format(response.status_code, response.reason))
    data = response.json()
    if isinstance(data, dict):
        messages = data.pop("messages", None)  
        if messages:
            print(json.dumps(messages, indent=4))
    try:
        response.raise_for_status()
    except Exception as e:
        print(e)
        return None
    else:
        return data

def search_stac(latlim, lonlim, timelim, product_name, extra_search_kwargs):

    timelim = adjust_timelim_dtype(timelim)
    sd = datetime.datetime.strftime(timelim[0], "%Y-%m-%dT00:00:00Z")
    ed = datetime.datetime.strftime(timelim[1], "%Y-%m-%dT23:59:59Z")
    search_dates = f"{sd}/{ed}"

    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    platform = {"LC08": "LANDSAT_8"}[product_name]
    search_kwargs = {   
                        **{'platform': {'or':[platform]}}, # TODO this doesnt work, so filtering manually later, doc example: `'platform': {'or':['LANDSAT_8','LANDSAT_9']}`.
                        **extra_search_kwargs
                    }

    stac = 'https://landsatlook.usgs.gov/stac-server' # Landsat STAC API Endpoint
    stac_response = r.get(stac).json() 
    catalog_links = stac_response['links']
    search = [l['href'] for l in catalog_links if l['rel'] == 'search'][0]   #retreive search endpoint from STAC Catalog

    params = dict()
    params['collections'] = ['landsat-c2l2-sr','landsat-c2l2-st'] 
    params['limit'] = 400
    params['bbox'] = bb
    params['datetime'] = search_dates
    params['query'] = search_kwargs    

    query = r.post(search, json=params, ).json()   # send POST request to the stac-search endpoint with params passed in

    ids = np.unique([x["id"].replace("_ST", "").replace("_SR", "") for x in query["features"] if product_name in x["id"]]).tolist()

    log.info(f"--> Found {len(ids)} scenes.")

    return ids, query

def request_scenes(ids, image_extents):

    uauth = pywapor.collect.accounts.get("EARTHEXPLORER") 

    order = espa_api('available-products', 
                        body = {"inputs": ids}, 
                        uauth = uauth)

    if "not implemented" in order.keys():
        missing = order.pop("not implemented")
        log.warning(f"--> Some scenes could not be found (`{'`, `'.join(missing)}`).")

    req_prods = ["l1"]
    for sensor in order.keys():
        order[sensor]["products"] = [x for x in req_prods if x in order[sensor]["products"]]

    order['projection'] = {"lonlat": None}
    order['format'] = 'netcdf'
    order['resampling_method'] = 'nn'
    order['note'] = f'pyWaPOR_{pywapor.__version__}'
    order["image_extents"] = image_extents

    print(f"--> Placing order for {len(ids)} scenes.")

    order_response = espa_api('order', verb='post', body = order, uauth = uauth)

    return order_response

def unpack(fp, folder):
    fn, _ = os.path.splitext(os.path.split(fp)[-1])
    subfolder = os.path.join(folder, fn)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)
    shutil.unpack_archive(fp, subfolder)

def download_scenes(scene_ids, product_folder, latlim, lonlim):

    image_extents = {
            "east": lonlim[1],
            "north": latlim[1],
            "south": latlim[0],
            "west": lonlim[0],
            "units": "dd",
        }

    uauth = pywapor.collect.accounts.get("EARTHEXPLORER") 

    # Make sure the total length is not 0 so the while-loop is started.
    to_request, to_wait, to_download = ([True], [True], [True])

    while len(to_download) + len(to_wait) + len(to_request) > 0:

        # Check which scenes already exist locally
        paths = glob.glob(os.path.join(product_folder, "**", "*.nc"), recursive = True)
        unpacked_scenes = set([os.path.splitext(os.path.split(f)[-1])[0] for f in paths if re.search(r'L.{3}_.{4}_\d{6}_\d{8}_\d{8}_\d{2}_T\d.nc', f)])
        packed_scenes = {[y.name for y in tarfile.open(x, encoding='utf-8').getmembers() if ".nc" in y.name][0][:-3]: x for x in glob.glob(os.path.join(product_folder, "**", "*.tar.gz"), recursive = True)}
        to_unpack = [fp for k, fp in packed_scenes.items() if (k in scene_ids) and (k not in unpacked_scenes)]
        for fp in to_unpack:
            unpack(fp, product_folder)
        available_scenes = unpacked_scenes.union(set(packed_scenes.keys()))

        if np.all([x in available_scenes for x in scene_ids]):
            break

        to_request = set()
        to_wait = set()
        to_download = dict()

        # Get an overview of the existing orders.
        all_orders = espa_api(f"item-status", uauth = uauth)
        for order_id, order in all_orders.items():

            # Get specific order details.
            order_details = espa_api(f"order/{order_id}", uauth = uauth)
            if order_details["product_opts"].get("image_extents") != image_extents:
                continue

            for scene in order:

                if (scene["status"] == "complete") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_download[scene["name"]] = scene["product_dload_url"]
                    to_request.discard(scene["name"])
                    to_wait.discard(scene["name"])
                    print(f"--> {scene['name']} status: complete.")
                elif (scene["status"] == "oncache") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: oncache.")
                elif (scene["status"] == "onorder") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: onorder.")
                elif (scene["status"] == "queued") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: queued.")
                elif (scene["status"] == "processing") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: processing.")
                elif (scene["status"] == "error") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_wait.discard(scene["name"])
                    print(f"--> {scene['name']} status: error ({scene['note']}).")
                elif (scene["status"] == "retry") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: retry.")
                elif (scene["status"] == "unavailable") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.discard(scene["name"])
                    to_wait.add(scene["name"])
                    print(f"--> {scene['name']} status: unavailable.")
                elif (scene["status"] == "cancelled") and (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.add(scene["name"])
                    to_wait.discard(scene["name"])
                    print(f"--> {scene['name']} status: cancelled.")
                elif (scene["name"] in scene_ids) and (scene["name"] not in available_scenes):
                    to_request.add(scene["name"])
                    to_wait.discard(scene["name"])
                    print(f"--> {scene['name']} status: {scene['status']}.")
                else:
                    ...

        # Request missing scenes on (ESPA)
        if len(to_request) > 0:
            _ = request_scenes(to_request, image_extents)
        
        # Download completed scenes (ESPA)
        if len(to_download) > 0:
            fps = download_urls(list(to_download.values()), product_folder)

            for fp in fps:
                unpack(fp, product_folder)
    
    return available_scenes

def _process_scene(scene, product_folder, variables, example_ds = None):

    ds = xr.open_dataset(scene, mask_and_scale=False, chunks = "auto")[list(variables.keys())]

    if len(set([ds[x].shape for x in ds.data_vars])) > 1:
        log.warning("--> Not all variables have identical shapes.")

    xdim = ds[list(variables.values())[0][0][1]].values
    ydim = ds[list(variables.values())[0][0][0]].values

    renames1 = {k: v[1] for k, v in variables.items()}
    renames2 = {v[0][0]: "y" for v in variables.values()}
    renames3 = {v[0][1]: "x" for v in variables.values()}

    ds = ds.rename_dims({**renames2, **renames3})
    ds = ds.drop(list(ds.coords.keys())).assign_coords({"x": xdim, "y": ydim})
    ds = ds.rename(renames1)

    crs = rasterio.crs.CRS.from_epsg(4326)
    ds = ds.rio.write_crs(crs)
    ds = ds.rio.write_grid_mapping("spatial_ref")

    ds = ds.sortby("y", ascending = False)
    ds = ds.sortby("x")

    mtl_fp = glob.glob(os.path.join(product_folder, "**", f"*{ds.LPGSMetadataFile.replace('.txt', '.xml')}*"), recursive=True)[0]
    with open(mtl_fp,"r") as f:
        xml_content = f.read()
    mtl = xmltodict.parse(xml_content)

    # Clip and pad to bounding-box
    if isinstance(example_ds, type(None)):
        example_ds = make_example_ds(ds, product_folder, crs, bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]])
    ds = ds.rio.reproject_match(example_ds).chunk("auto")
    ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

    # Apply variable specific functions.
    for vars in variables.values():
        for func in vars[2]:
            ds, label = apply_enhancer(ds, vars[1], func)
            print(label)

    # Add time dimension to data arrays.
    ds = ds.expand_dims({"time": 1})

    # Set the correct time.
    date_str = mtl['LANDSAT_METADATA_FILE']["IMAGE_ATTRIBUTES"]["DATE_ACQUIRED"]
    time_str = mtl['LANDSAT_METADATA_FILE']["IMAGE_ATTRIBUTES"]["SCENE_CENTER_TIME"]
    datetime_str = date_str + " " + time_str.replace("Z", "")
    ds = ds.assign_coords({"time":[np.datetime64(datetime_str)]})

    # Cleanup attributes
    for var in ["x", "y", "time"]:
        ds[var].attrs = {}

    return ds, example_ds

def process_scenes(fp, scene_paths, product_folder, variables, post_processors):
    
    dss = list()

    log.info(f"--> Processing {len(scene_paths)} scenes.").add()

    example_ds = None
    for i, scene in enumerate(scene_paths):
        ds, example_ds = _process_scene(scene, product_folder, variables, example_ds = example_ds)
        fp_temp = scene.replace(".nc", "_temp.nc")
        ds = save_ds(ds, fp_temp, encoding = "initiate", label = f"({i+1}/{len(scene_paths)}) Processing `{os.path.split(scene)[-1]}`.")
        dss.append(ds)

    ds = xr.concat(dss, "time")
    
    # Apply general product functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    # Remove unrequested variables.
    ds = ds[list(post_processors.keys())]
    
    for var in ds.data_vars:
        ds[var].attrs = {}

    ds = ds.sortby("time")

    # Save final netcdf.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
        warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
        ds = save_ds(ds, fp, chunks = "auto", encoding = "initiate", label = f"Merging files.")

    # Remove intermediate files.
    for x in dss:
        remove_ds(x)

    return ds

def download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = None, post_processors = None, 
                extra_search_kwargs = {'eo:cloud_cover': {'gte': 0, 'lt': 30}}):

    product_folder = os.path.join(folder, "LANDSAT")

    appending = False
    fn = os.path.join(product_folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        os.rename(fn, fn.replace(".nc", "_to_be_appended.nc"))
        existing_ds = open_ds(fn.replace(".nc", "_to_be_appended.nc"))
        if np.all([x in existing_ds.data_vars for x in req_vars]):
            existing_ds = existing_ds.close()
            os.rename(fn.replace(".nc", "_to_be_appended.nc"), fn)
            existing_ds = open_ds(fn)
            return existing_ds[req_vars]
        else:
            appending = True
            fn = os.path.join(product_folder, f"{product_name}_appendix.nc")
            req_vars = [x for x in req_vars if x not in existing_ds.data_vars]

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items() if k in req_vars}

    # Search scene IDs (STAC)
    scene_ids, _ = search_stac(latlim, lonlim, timelim, product_name, extra_search_kwargs)

    # Order and download scenes (ESPA)
    available_scenes = download_scenes(scene_ids, product_folder, latlim, lonlim)

    # Process scenes.
    scene_paths = [glob.glob(os.path.join(product_folder, "**", f"*{x}.nc"), recursive = True)[0] for x in available_scenes]
    ds_new = process_scenes(fn, scene_paths, product_folder, variables, post_processors)

    if appending:
        ds = xr.merge([ds_new, existing_ds])
        lbl = f"Appending new variables (`{'`, `'.join(req_vars)}`) to existing file."
        ds = save_ds(ds, os.path.join(product_folder, f"{product_name}.nc"), encoding = "initiate", label = lbl)
        remove_ds(ds_new)
        remove_ds(existing_ds)
    else:
        ds = ds_new

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Local/landsat_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2022-03-29", "2022-04-25"]
    product_name = "LC08"        
    req_vars = ["red", "nir", "lst"]
    variables = None
    post_processors = None
    example_ds = None
    extra_search_kwargs = {'eo:cloud_cover': {'gte': 0, 'lt': 30}}

    scene = r"/Users/hmcoerver/Local/landsat_test/LANDSAT/LC081770402022041302T1-SC20221214143308.tar/LC08_L2SP_177040_20220413_20220420_02_T1.nc"
    product_folder = os.path.join(folder, "LANDSAT")
    variables = default_vars(product_name, req_vars)
