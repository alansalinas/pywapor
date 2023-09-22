"""
Code to download data using the Copernicus OData API.
https://documentation.dataspace.copernicus.eu/APIs/OData.html
"""
import time
import re
import tqdm
import glob
import requests
import datetime
import os
import shutil
import warnings
import rasterio.crs
import numpy as np
import xarray as xr
from cachetools import cached, FIFOCache, TTLCache
from pywapor.general.logger import adjust_logger, log
from pywapor.collect.protocol.crawler import download_urls
from pywapor.general.processing_functions import save_ds, open_ds, unpack, make_example_ds, remove_ds
from pywapor.general.logger import log
from pywapor.collect.accounts import get
from pywapor.enhancers.apply_enhancers import apply_enhancer

@cached(cache=TTLCache(maxsize=2048, ttl=550))
def get_access_token():
    """
    Create account at https://dataspace.copernicus.eu/
    """

    username, password = get("COPERNICUS_DATA_SPACE")

    data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password",
        }
    
    try:
        log.info("--> Requesting access token.")
        r = requests.post("https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
                            data=data,
                            )
        r.raise_for_status()
    except Exception as e:
        raise Exception(
                        f"Access token creation failed. Reponse from the server was: {r.json()}"
                    )
    
    token = r.json()

    return token

@cached(cache=FIFOCache(maxsize=2048))
def get_scene_nodes(scene_id):
    files = list()
    nodes_to_crawl = [f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products({scene_id})/Nodes"]
    
    while len(nodes_to_crawl) > 0:
        node_path = nodes_to_crawl[0]

        node_ = requests.get(node_path)
        node_.raise_for_status()
        node = node_.json()["result"]

        for sub_node in node:
            sub_node_name = sub_node["Name"]
            if sub_node["ChildrenNumber"] != 0:
                nodes_to_crawl.append(node_path + f"({sub_node_name})/Nodes")
            else:
                uri = sub_node["Nodes"]["uri"]
                files.append(uri[:-6] if uri[-6:] == "/Nodes" else uri)

        nodes_to_crawl.remove(node_path)
    
    return files

def download(folder, latlim, lonlim, timelim, product_name, product_type, node_filter = None):

    sd = datetime.datetime.strftime(timelim[0], "%Y-%m-%dT00:00:00Z")
    ed = datetime.datetime.strftime(timelim[1], "%Y-%m-%dT23:59:59Z")

    product_name_ = {"SENTINEL3": "SENTINEL-3", 
                    "SENTINEL2": "SENTINEL-2",
                    "SENTINEL-2": "SENTINEL-2",
                    "SENTINEL-3": "SENTINEL-3",
                    }[product_name]

    base_url = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
    filters = [
                f"Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and att/OData.CSC.StringAttribute/Value eq '{product_type}')",
                f"OData.CSC.Intersects(area=geography'SRID=4326;POLYGON(({lonlim[0]} {latlim[0]},{lonlim[1]} {latlim[0]},{lonlim[1]} {latlim[1]},{lonlim[0]} {latlim[1]},{lonlim[0]} {latlim[0]}))')",
                f"Collection/Name eq '{product_name_}'",
                f"ContentDate/Start gt {sd}",
                f"ContentDate/Start lt {ed}",
            ]

    query = base_url + "?$filter=" + " and ".join(filters)

    results = {"@odata.nextLink": query}
    scenes = list()    
    while "@odata.nextLink" in results.keys():
        out = requests.get(results["@odata.nextLink"])
        out.raise_for_status()
        results = out.json()
        scenes += results["value"]

    # Drop identical scenes.
    scene_names = {x["Name"]: i for i, x in enumerate(scenes)}
    scenes = [scenes[i] for i in scene_names.values()]

    # Filter reprocessed scenes (select newest).
    overview = dict()
    for i, scene in enumerate(scenes):
        mission_id, level, stime, baseline, rel_orbit, tile_number, rest = scene["Name"].split("_")
        key = "_".join([mission_id, level, stime, rel_orbit, tile_number])
        baseline_int = int(baseline[1:])
        if not key in list(overview.keys()):
            overview[key] = (baseline_int, scene)
        else:
            if baseline_int > overview[key][0]:
                overview[key] = (baseline_int, scene)

    scenes = [x[1] for x in overview.values()]

    log.info(f"--> Searching nodes for {len(scenes)} `{product_name_}.{product_type}` scenes.")

    urls = list()
    fps = list()
    fp_checks = list()

    for scene in tqdm.tqdm(scenes, leave = False):
        nodes = get_scene_nodes(scene["Id"])
        if not isinstance(node_filter, type(None)):
            filtered_nodes = set([x for x in nodes if node_filter(x)])
        else:
            filtered_nodes = set(nodes)

        for node in filtered_nodes:
            node_path_rel = re.findall("Nodes\((.*?)\)", node)
            fp = os.path.join(folder, *node_path_rel)
            url = node + "/$value"
            fp_checks.append(fp)
            if not os.path.isfile(fp):
                urls.append(url)
                fps.append(fp)

    log.info(f"--> {len(fp_checks)} nodes required.")
    log.info(f"--> Downloading {len(urls)} missing nodes.").add()

    block_size = 10
    dled_fps = list()
    for i in range(0, len(urls), block_size):
        token = get_access_token()
        access_token = token["access_token"]
        headers = {'Authorization': f'Bearer {access_token}'}
        fps_ = download_urls(urls[i:i + block_size], None, fps = fps[i:i + block_size], headers = headers)
        dled_fps += fps_
        time.sleep(30)

    log.sub()

    for fp_ in fp_checks:
        if not os.path.isfile(fp_):
            log.warning(f"--> `{fp}` missing.")

    if "3" in product_name: # NOTE this is lazy...
        dled_scenes = sorted(list(set([re.compile(".*\.SEN3").search(x).group() for x in fp_checks])))
    if "2" in product_name:
        dled_scenes = sorted(list(set([re.compile(".*\.SAFE").search(x).group() for x in fp_checks])))

    return list(dled_scenes)

def process_sentinel(scenes, variables, time_func, final_fn, post_processors, processor, bb = None):
    """Process downloaded Sentinel scenes into netCDFs.

    Parameters
    ----------
    scenes : list
        Paths to downloaded nodes.
    variables : dict
        Keys are variable names, values are additional settings.
    time_func : function
        Function that parses a np.datetime64 from a filename.
    final_fn : str
        Path to the file in which to store all the combined data.
    post_processors : dict
        Functions to apply when the data has been processed.
    processor : function
        Function to apply sensor specific transformations.
    bb : list, optional
        Boundingbox to clip to, [xmin, ymin, xmax, ymax], by default None.

    Returns
    -------
    xr.Dataset
        Ouput data.

    Raises
    ------
    ValueError
        Invalid value for `source_name`.
    """

    chunks = {"time": 1, "x": 1000, "y": 1000}

    example_ds = None
    dss1 = dict()

    log.info(f"--> Processing {len(scenes)} scenes.").add()

    target_crs = rasterio.crs.CRS.from_epsg(4326)

    for i, scene_folder in enumerate(scenes):
        
        folder, fn = os.path.split(scene_folder)

        ext = os.path.splitext(scene_folder)[-1]
        
        fp = os.path.join(folder, os.path.splitext(fn)[0] + ".nc")
        if os.path.isfile(fp):
            ds = open_ds(fp)
            dtime = ds.time.values[0]
            if dtime in dss1.keys():
                dss1[dtime].append(ds)
            else:
                dss1[dtime] = [ds]
            continue

        if ext == ".zip":
            scene_folder = unpack(fn, folder)
            remove_folder = True
        else:
            scene_folder = scene_folder
            remove_folder = False

        ds = processor(scene_folder, variables, bb = bb)

        # Apply variable specific functions.
        for vars in variables.values():
            for func in vars[2]:
                ds, label = apply_enhancer(ds, vars[1], func)

        # NOTE: see https://github.com/corteva/rioxarray/issues/545
        # NOTE: see rioxarray issue here: https://github.com/corteva/rioxarray/issues/570
        ds = ds.sortby("y", ascending = False)
        _ = [ds[var].rio.write_nodata(np.nan, inplace = True) for var in ds.data_vars]

        # Clip and pad to bounding-box
        if isinstance(example_ds, type(None)):
            example_ds = make_example_ds(ds, folder, target_crs, bb = bb)
        ds = ds.rio.reproject_match(example_ds).chunk({k: v for k, v in chunks.items() if k in ["x", "y"]})

        # NOTE: see rioxarray issue here: https://github.com/corteva/rioxarray/issues/570
        _ = [ds[var].attrs.pop("_FillValue") for var in ds.data_vars if "_FillValue" in ds[var].attrs.keys()]
        
        ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

        dtime = time_func(fn)
        ds = ds.expand_dims({"time": 1})
        ds = ds.assign_coords({"time":[dtime]})
        ds.attrs = {}

        # Save to netcdf
        ds = save_ds(ds, fp, chunks = chunks, encoding = "initiate", label = f"({i+1}/{len(scenes)}) Processing {fn} to netCDF.")

        if dtime in dss1.keys():
            dss1[dtime].append(ds)
        else:
            dss1[dtime] = [ds]

        if remove_folder:
            shutil.rmtree(scene_folder)

    log.sub()
    
    # Merge spatially.
    dss = [xr.concat(dss0, "stacked").median("stacked") for dss0 in dss1.values()]

    # Merge temporally.
    ds = xr.concat(dss, "time")

    # Define output path.
    fp = os.path.join(folder, final_fn)
    
    # Apply general product functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    # Remove unrequested variables.
    ds = ds[list(post_processors.keys())]
    
    for var in ds.data_vars:
        ds[var].attrs = {}

    ds = ds.rio.write_crs(target_crs)

    ds = ds.sortby("time")

    # Save final netcdf.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
        warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
        ds = save_ds(ds, fp, chunks = chunks, encoding = "initiate", label = f"Merging files.")

    # Remove intermediate files.
    for dss0 in dss1.values():
        for x in dss0:
            remove_ds(x)

    return ds

if __name__ == "__main__":

    product_name = "SENTINEL-2"
    product_type = "S2MSI2A"

    # product_name = "SENTINEL-3"
    # product_type = "SL_2_LST___"

    timelim = [datetime.date(2023, 3, 1), datetime.date(2023, 3, 3)]
    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]
    folder = r"/Users/hmcoerver/Local/new_sentinel_test"

    adjust_logger(True, folder, "INFO")

    variables = {
        "_B01_60m.jp2": [(), "coastal_aerosol", []],
        "_B02_60m.jp2": [(), "blue", []],
        "_B03_60m.jp2": [(), "green", []],
        "_B04_60m.jp2": [(), "red", []],
    }

    # variables = {
    #     "LST_in.nc": [(), "lst", []],
    #     "geodetic_in.nc": [(), "coords", []],
    # }

    def node_filter(node_info):
        fn = os.path.split(node_info)[-1]
        to_dl = list(variables.keys()) + ["MTD_MSIL2A.xml"]
        return np.any([x in fn for x in to_dl])
    
    dled_scenes = download(folder, latlim, lonlim, timelim, product_name, product_type, node_filter = node_filter)

