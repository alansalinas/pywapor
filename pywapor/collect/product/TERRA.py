import datetime
import requests
import os
import copy
import numpy as np
import xarray as xr
import pywapor.collect.protocol.cog as cog
import pywapor.general.bitmasks as bm
import pywapor.collect.accounts as accounts
from functools import partial
from osgeo import gdal
from pywapor.general.logger import log, adjust_logger
from pywapor.general.processing_functions import save_ds, open_ds, remove_ds, adjust_timelim_dtype

def default_post_processors(product_name, req_vars = ["ndvi", "r0"]):
    """Given a `product_name` and a list of requested variables, returns a dictionary with a 
    list of functions per variable that should be applied after having collected the data
    from a server.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list, optional
        List of variables to be collected, by default ["ndvi", "r0"].

    Returns
    -------
    dict
        Functions per variable that should be applied to the variable.
    """

    post_processors = {
        "urn:eop:VITO:PROBAV_S5_TOC_100M_COG_V2": {
            "r0": [
                    calc_r0, 
                    partial(mask_bitwise_qa, 
                    flags = ["bad BLUE", "bad RED", "bad NIR", "bad SWIR", 
                            "sea", "undefined", "cloud", "ice/snow", "shadow"]),
            ],
            "ndvi": [
                    partial(mask_bitwise_qa, 
                    flags = ["bad RED", "bad NIR", "sea", "undefined", 
                            "cloud", "ice/snow", "shadow"]),
            ],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def default_vars(product_name, req_vars = ["ndvi", "r0"]):
    """Given a `product_name` and a list of requested variables, returns a dictionary
    with metadata on which exact layers need to be requested from the server, how they should
    be renamed, and how their dimensions are defined.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list, optional
        List of variables to be collected, by default ["ndvi", "r0"].

    Returns
    -------
    dict
        Metadata on which exact layers need to be requested from the server.
    """
    
    variables = {
        "urn:eop:VITO:PROBAV_S5_TOC_100M_COG_V2": {
            "NIR":      [("lon", "lat"), "nir"],
            "RED":      [("lon", "lat"), "red"],
            "BLUE":     [("lon", "lat"), "blue"],
            "SWIR":     [("lon", "lat"), "swir"],
            "SM":       [("lon", "lat"), "qa"],
            "GEOMETRY": [("lon", "lat"), "geometry"],
            "NDVI":     [("lon", "lat"), "ndvi"],
            "TIME":     [("lon", "lat"), "time"],
        },
    }

    req_dl_vars = {
        "urn:eop:VITO:PROBAV_S5_TOC_100M_COG_V2": {
            "ndvi": ["NDVI", "SM"],
            "r0":   ["BLUE", "NIR", "RED", "SWIR", "SM"],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def mask_bitwise_qa(ds, var, flags):
    """Mask PROBAV data using a qa variable.

    Parameters
    ----------
    ds : xr.Dataset
        Input data.
    var : str
        Variable in `ds` to mask.
    flags : list
        Which flags not to mask.

    Returns
    -------
    xr.Dataset
        Masked dataset.
    """
    flag_bits = bm.PROBAV_qa_translator()
    mask = bm.get_mask(ds["qa"].astype("uint8"), flags, flag_bits)
    ds[var] = ds[var].where(~mask, np.nan)
    return ds

def calc_r0(ds, *args):
    ds["r0"] = 0.429 * ds["blue"] + 0.333 * ds["red"] + 0.133 * ds["nir"] + 0.105 * ds["swir"]
    return ds

def create_header(folder):

    un, pw = accounts.get("TERRA")

    url = "https://sso.vgt.vito.be/auth/realms/terrascope/protocol/openid-connect/token"
    params = {
        "grant_type": "password",
        "client_id": "public",
        "username": un,
        "password": pw,
    }
    headers = {'content-type': 'application/x-www-form-urlencoded'}

    resp = requests.post(url, data = params, headers = headers)
    resp.raise_for_status()
    token = resp.json()["access_token"]

    header_file = os.path.join(folder, "token.txt")
    if os.path.isfile(header_file):
        try:
            os.remove(header_file)
        except PermissionError:
            log.warning("--> Unable to remove `token.txt` file (PermissionError).")


    with open(header_file, "w") as text_file:
        text_file.write(f"Authorization: Bearer {token}")

    return header_file

def download(folder, latlim, lonlim, timelim, product_name, req_vars = ["ndvi", "r0"],
                variables = None, post_processors = None, timedelta = np.timedelta64(60, "h")):
    """Download MODIS data and store it in a single netCDF file.

    Parameters
    ----------
    folder : str
        Path to folder in which to store results.
    latlim : list
        Latitude limits of area of interest.
    lonlim : list
        Longitude limits of area of interest.
    timelim : list
        Period for which to prepare data.
    product_name : str
        Name of the product to download.
    req_vars : list, optional
        Which variables to download for the selected product, by default ["ndvi", "r0"].
    variables : dict, optional
        Metadata on which exact layers need to be requested from the server, by default None.
    post_processors : dict, optional
        Functions per variable that should be applied to the variable, by default None.

    Returns
    -------
    xr.Dataset
        Downloaded data.
    """
    folder = os.path.join(folder, "TERRA")

    if not os.path.isdir(folder):
        os.makedirs(folder)

    fn = os.path.join(folder, f"{product_name.replace(':','_')}.nc")
    req_vars_orig = copy.deepcopy(req_vars)
    if os.path.isfile(fn):
        existing_ds = open_ds(fn)
        req_vars_new = list(set(req_vars).difference(set(existing_ds.data_vars)))
        if len(req_vars_new) > 0:
            req_vars = req_vars_new
            existing_ds = existing_ds.close()
        else:
            return existing_ds[req_vars_orig]
        
    spatial_buffer = True
    if spatial_buffer:
        dx = dy = 0.0033
        latlim = [latlim[0] - dy, latlim[1] + dy]
        lonlim = [lonlim[0] - dx, lonlim[1] + dx]

    timelim = adjust_timelim_dtype(timelim)

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items() if k in req_vars}

    bb = [lonlim[0],latlim[0],lonlim[1],latlim[1]]

    coords = {"x": ["lon", lonlim], "y": ["lat", latlim]}

    sd = datetime.datetime.strftime(timelim[0], "%Y-%m-%dT00:00:00Z")
    ed = datetime.datetime.strftime(timelim[1], "%Y-%m-%dT23:59:59Z")
    search_dates = f"{sd}/{ed}"
    params = dict()
    params['limit'] = 250
    params['bbox'] = bb
    params['datetime'] = search_dates
    params['collections'] = [product_name]

    products = search_stac(params)

    dss = list()

    log.info(f"--> Processing {len(products)} scenes.").add()

    for i, product in enumerate(products):

        out_fp = os.path.join(folder, f"{product['properties']['title']}.nc")

        header_file = create_header(folder)
        all_cog_urls = [v.get("href") for k, v in product["assets"].items()]
        cog_urls = [all_cog_urls[np.argmax([key in x for x in all_cog_urls])] for key in variables.keys()]
        dl_urls = [f"/vsicurl?header_file={header_file}&url={cog_url}" for cog_url in cog_urls]

        ds_ = gdal.BuildVRT(out_fp.replace(".nc", ".vrt"), dl_urls, options = gdal.BuildVRTOptions(separate = True))
        ds_.FlushCache()
        ds_ = None

        variables_ = {f"Band{np.argmax([key in x for x in dl_urls]) + 1}": variables[key] for key in variables.keys()}

        log.info(f"--> ({i+1}/{len(products)}) Downloading `{product['properties']['title']}`.").add()

        ds = cog.download(out_fp, product_name, coords, variables_, post_processors=post_processors, url_func=lambda x: out_fp.replace(".nc", ".vrt"))
        ds = ds.expand_dims({"time": 1}).assign_coords({"time": [datetime.datetime.strptime(product["properties"]["datetime"], "%Y-%m-%dT%H:%M:%SZ")]})

        dss.append(ds)

        log.sub()
    
    log.sub()

    ds = xr.merge(dss)

    if not isinstance(timedelta, type(None)):
        ds["time"] = ds["time"] + timedelta

    ds = ds[req_vars]
    ds = ds.rio.write_crs("epsg:4326")

    ds = save_ds(ds, fn, label = f"Merging files.")

    for nc in dss:
        remove_ds(nc)

    return ds[req_vars_orig]

def search_stac(params):
    endpoint = 'https://services.terrascope.be/stac'
    stac_response = requests.get(endpoint)
    stac_response.raise_for_status()
    catalog_links = stac_response.json()['links']
    search = [l['href'] for l in catalog_links if l['rel'] == 'search'][0]

    all_scenes = list()

    links = [{"body": params, "rel": "next"}]

    while "next" in [x.get("rel", None) for x in links]:
        params_ = [x.get("body") for x in links if x.get("rel") == "next"][0]
        query = requests.post(search, json = params_)
        query.raise_for_status()
        out = query.json()
        all_scenes += out["features"]
        links = out["links"]

    return all_scenes
    
if __name__ == "__main__":

    product_name = "urn:eop:VITO:PROBAV_S5_TOC_100M_COG_V2"

    folder = r"/Users/hmcoerver/Local/bla_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 2), datetime.date(2020, 7, 12)]

    adjust_logger(True, folder, "INFO")

    variables = None
    post_processors = None
    req_vars = ["r0", "ndvi"]