import requests
import datetime
import base64
import copy
import os
import xarray as xr
import numpy as np
from osgeo import gdal
from functools import partial
from pywapor.collect import accounts
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.processing_functions import save_ds, remove_ds, open_ds
from pywapor.general.logger import log, adjust_logger
from pywapor.enhancers.other import drop_empty_times
from pywapor.general.curvilinear import curvi_to_recto
gdal.UseExceptions()

def get_token():
    """_summary_

    Returns
    -------
    _type_
        _description_
    """
    s3_endpoint = "https://data.laadsdaac.earthdatacloud.nasa.gov/s3credentials"
    
    login_resp = requests.get(s3_endpoint, allow_redirects=False)
    login_resp.raise_for_status()

    un, pw = accounts.get("NASA")
    auth = f"{un}:{pw}"
    encoded_auth  = base64.b64encode(auth.encode('ascii'))

    auth_redirect = requests.post(
        login_resp.headers['location'],
        data = {"credentials": encoded_auth},
        headers= {"Origin": s3_endpoint},
        allow_redirects=False
    )
    auth_redirect.raise_for_status()
    if "currently unavailable" in auth_redirect.text:
        log.warning(f"The Earthdata Service is currently unavailable. Check {login_resp.headers['location']}.")
        raise ConnectionError(f"Earthdata Service is currently unavailable. Check {login_resp.headers['location']}.")

    final = requests.get(auth_redirect.headers['location'], allow_redirects=False)
    results = requests.get(s3_endpoint, cookies={'accessToken': final.cookies['accessToken']})
    results.raise_for_status()

    return results.json()
    
def download_arrays(path, subdss, folder, path_appendix = ".nc", creation_options = ["ARRAY:COMPRESS=DEFLATE"]):

    fn = os.path.splitext(os.path.split(path)[-1])[0]
    path_local = os.path.join(folder, fn + path_appendix)

    if os.path.isfile(path_local):
        ...
    else:
        creds = get_token()

        gdal_config_options = {
            "AWS_ACCESS_KEY_ID": creds["accessKeyId"],
            "AWS_SESSION_TOKEN": creds["sessionToken"],
            "AWS_SECRET_ACCESS_KEY": creds["secretAccessKey"],
            "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        }

        try:
            for k, v in gdal_config_options.items():
                gdal.SetConfigOption(k, v)

            md_options = gdal.MultiDimTranslateOptions(
                    format = "netCDF",
                    arraySpecs = subdss,
                    creationOptions = creation_options
            )

            ds = gdal.MultiDimTranslate(path_local, path, options = md_options)
            ds = ds.FlushCache()
        except Exception as e:
            raise e
        finally:
            for k, v in gdal_config_options.items():
                gdal.SetConfigOption(k, None)
        
    return path_local

def search_stac(params, endpoint = 'https://cmr.earthdata.nasa.gov/stac/LAADS'):
    stac_response = requests.get(endpoint)
    stac_response.raise_for_status()
    catalog_links = stac_response.json()['links']
    search = [l['href'] for l in catalog_links if l['rel'] == 'search'][0]
    query = requests.post(search, json=params)
    query.raise_for_status()
    out = query.json()
    if out["context"]["returned"] < out["context"]["matched"]:
        log.warning("Number of matched features exceeds limit.") # TODO: make warning
    return out["features"]

def create_stac_summary(bb, timelim):
    
    sd = datetime.datetime.strftime(timelim[0], "%Y-%m-%dT00:00:00Z")
    ed = datetime.datetime.strftime(timelim[1], "%Y-%m-%dT23:59:59Z")
    search_dates = f"{sd}/{ed}"

    params = dict()
    params['limit'] = 250
    params['bbox'] = bb
    params['datetime'] = search_dates

    params['collections'] = ['VNP02IMG.v2']
    results02 = search_stac(params)
    params['collections'] = ['VNP03IMG.v2']
    results03 = search_stac(params)
    params['collections'] = ['CLDMSK_L2_VIIRS_SNPP.v1']
    resultsqa = search_stac(params)

    get_fns = lambda results: [os.path.split(x["assets"]["data"]["href"])[-1] for x in results]
    get_dts = lambda files: {datetime.datetime.strptime(".".join(x.split(".")[1:3]), "A%Y%j.%H%M"): x for x in files}
    files02 = get_fns(results02)
    files03 = get_fns(results03)
    filesqa = get_fns(resultsqa)

    dates02 = get_dts(files02)
    dates03 = get_dts(files03)
    datesqa = get_dts(filesqa)

    unique_dates = set(dates02.keys()).union(set(dates03.keys())).union(set(datesqa.keys()))

    summary = {k: (dates02.get(k, None), dates03.get(k, None), datesqa.get(k, None)) for k in unique_dates}
    filtered_summary = {k: v for k, v in summary.items() if not None in  v}

    return filtered_summary

def default_vars(product_name, req_vars):
    """Given a `product_name` and a list of requested variables, returns a dictionary
    with metadata on which exact layers need to be requested from the server, how they should
    be renamed, and how their dimensions are defined.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list
        List of variables to be collected.

    Returns
    -------
    dict
        Metadata on which exact layers need to be requested from the server.
    """

    variables = {
        "VNP02IMG": {"bt": {(), "bt"}},
    }

    req_dl_vars = {
        "VNP02IMG": {
            "bt": ["bt"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars = None):
    """Given a `product_name` and a list of requested variables, returns a dictionary with a 
    list of functions per variable that should be applied after having collected the data
    from a server.

    Parameters
    ----------
    product_name : str
        Name of the product.
    req_vars : list
        List of variables to be collected.

    Returns
    -------
    dict
        Functions per variable that should be applied to the variable.
    """
    post_processors = {
        "VNP02IMG": {
            "bt": [partial(drop_empty_times, drop_vars = ["bt"])]
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def preproc(ds):
    ds = ds.rename({"lat": "y", "lon": "x", "Band1": "bt"})
    ds = ds.drop("crs")
    ds.attrs = {}
    for var in ds.data_vars:
        ds[var].attrs = {}
    ds = ds.rio.write_grid_mapping("spatial_ref")
    ds = ds.rio.write_crs("epsg:4326")
    ds = ds.rio.write_transform(ds.rio.transform(recalc=True))
    fn = os.path.split(ds.encoding["source"])[-1]
    date = np.datetime64(fn[3:13] + " " + fn[14:22].replace("_", ":"))
    return ds.expand_dims("time").assign_coords({"time": [date]})

def combine_unprojected_data(nc02_file, ncqa_file, lut_file, unproj_fn):

    ds1 = xr.open_dataset(nc02_file, mask_and_scale=False)
    ds2 = xr.open_dataset(lut_file, mask_and_scale=False)
    ds3 = xr.open_dataset(ncqa_file, chunks = "auto")

    # Create cloud mask. (0=cloudy, 1=probably cloudy, 2=probably clear, 3=confident clear, -1=no result)
    cmask = ds3["Integer_Cloud_Mask"].interp({
                    "number_of_lines": np.linspace(0-0.25, (ds3["number_of_lines"].size-1)+0.25, ds3["number_of_lines"].size*2),
                    "number_of_pixels": np.linspace(0-0.25, (ds3["number_of_pixels"].size-1)+0.25, ds3["number_of_pixels"].size*2)
                    }, kwargs = {"fill_value": "extrapolate"}).drop_vars(["number_of_lines", "number_of_pixels"])

    # Convert DN to Brightness Temperature using the provided lookup table.
    bt_da = ds2["I05_brightness_temperature_lut"].isel(number_of_LUT_values = ds1["I05"])

    # Chunk and mask invalid pixels.
    bt = bt_da.chunk("auto").where((bt_da >= bt_da.attrs["valid_min"]) & 
                                        (bt_da <= bt_da.attrs["valid_max"]) & 
                                        (bt_da != bt_da.attrs["_FillValue"]) &
                                        (ds1["I05_quality_flags"] == 0) &
                                        (cmask == 3), bt_da.attrs["_FillValue"]
    ).rename("bt").to_dataset()

    _ = save_ds(bt, unproj_fn, encoding = "initiate", label = "Combining data.").close()

def create_grid(latlim, lonlim, dx_dy = (0.0033, 0.0033)):
    dx, dy = dx_dy
    nx = np.ceil((lonlim[1] - lonlim[0]) / dx)
    ny = np.ceil((latlim[1] - latlim[0]) / dy)
    bb = [lonlim[0], latlim[0], lonlim[0] + nx * dx, latlim[0] + ny * dy]
    return bb, nx, ny

def download(folder, latlim, lonlim, timelim, product_name, req_vars,
                variables = None, post_processors = None):
    """Download VIIRSL1 data and store it in a single netCDF file.

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
    req_vars : list
        Which variables to download for the selected product.
    variables : dict, optional
        Metadata on which exact layers need to be requested from the server, by default None.
    post_processors : dict, optional
        Functions per variable that should be applied to the variable, by default None.

    Returns
    -------
    xr.Dataset
        Downloaded data.
    """
    folder = os.path.join(folder, "VIIRSL1")
    if not os.path.exists(folder):
        os.makedirs(folder)

    fn = os.path.join(folder, f"{product_name}.nc")
    req_vars_orig = copy.deepcopy(req_vars)
    if os.path.isfile(fn):
        existing_ds = open_ds(fn)
        req_vars_new = list(set(req_vars).difference(set(existing_ds.data_vars)))
        if len(req_vars_new) > 0:
            req_vars = req_vars_new
            existing_ds = existing_ds.close()
        else:
            return existing_ds[req_vars_orig]

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items() if k in req_vars}

    bb, nx, ny = create_grid(latlim, lonlim, dx_dy = (0.0033, 0.0033))
    filtered_summary = create_stac_summary(bb, timelim)

    log.info(f"--> Found {len(filtered_summary)} VIIRS scenes.").add()

    buckets = ("VNP02IMG", "VNP03IMG","CLDMSK_L2_VIIRS_SNPP")
    all_urls = {k: tuple(f"/vsis3/prod-lads/{buckets[i]}/{v[i]}" for i in range(3)) for k, v in filtered_summary.items()}
    all_proj_files = list()

    for i, (date, (nc02, nc03, nc_cloud)) in enumerate(all_urls.items()):

        date_str = str(date).replace(" ", "_").replace(":","_")
        unproj_fn = os.path.join(folder, f"bt_{date_str}_unprojected.nc")
        proj_fn = os.path.join(folder, f"bt_{date_str}_projected.nc")

        log.info(f"--> ({i+1}/{len(all_urls)}) Processing '{os.path.split(nc02)[-1]}'.").add()

        if os.path.isfile(proj_fn):
            all_proj_files.append(proj_fn)
            log.sub()
            continue

        lats_file = download_arrays(nc03, ["/geolocation_data/latitude"], folder, path_appendix="_lat.nc")
        lons_file = download_arrays(nc03, ["/geolocation_data/longitude"], folder, path_appendix="_lon.nc")
        nc02_file = download_arrays(nc02, ["/observation_data/I05_quality_flags", "/observation_data/I05"], folder)
        lut_file = download_arrays(nc02, ["/observation_data/I05_brightness_temperature_lut"], folder, path_appendix = "_lut.nc")
        ncqa_file = download_arrays(nc_cloud, ["/geophysical_data/Integer_Cloud_Mask"], folder)

        if not os.path.isfile(unproj_fn):
            combine_unprojected_data(nc02_file, ncqa_file, lut_file, unproj_fn)

        warp_kwargs = {"outputBounds": bb, "width": nx, "height": ny}
        _ = curvi_to_recto(lats_file, lons_file, {"bt": unproj_fn}, proj_fn, warp_kwargs = warp_kwargs)
        all_proj_files.append(proj_fn)

        for x in [nc02_file, ncqa_file, lut_file, unproj_fn]:
            remove_ds(x)

        log.sub()

    log.sub()

    ds = xr.open_mfdataset(all_proj_files, preprocess = preproc)

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    ds = save_ds(ds, fn, encoding = "initiate", label = "Merging files.")

    return ds[req_vars_orig]

if __name__ == "__main__":

    timelim = [datetime.datetime(2021, 3, 4), datetime.datetime(2021, 3, 5)]
    lonlim = [22.78125, 32.48438]
    latlim = [23.72777, 32.1612]
    folder = workdir = r"/Users/hmcoerver/Local/viirs"
    product_name = "VNP02IMG"
    req_vars = ["bt"]
    variables = None
    post_processors = None

    adjust_logger(True, folder, "INFO")

    ds = download(folder, latlim, lonlim, timelim, product_name, req_vars,
                    variables = variables, post_processors = post_processors)
        
        
