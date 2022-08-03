import xarray as xr
import numpy as np
import os
import tqdm
import rasterio
import rioxarray.merge
import tempfile
from pydap.cas.urs import setup_session
import urllib.parse
from pywapor.general.logger import log
from rasterio.crs import CRS
from pywapor.general.processing_functions import save_ds, process_ds
import warnings
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.collect.protocol.crawler import download_url, download_urls

def download(folder, product_name, coords, variables, post_processors, 
                fn_func, url_func, un_pw = None, tiles = None,  
                data_source_crs = None, parallel = False, spatial_tiles = True, 
                request_dims = True, timedelta = None):

    # Create selection object.
    selection = create_selection(coords, target_crs = data_source_crs)

    # Make output filepaths, should be same length as `urls`.
    fps = [os.path.join(folder, fn_func(product_name, x)) for x in tiles]

    # Make data request URLs.
    session = start_session(url_func(product_name, tiles[0]), selection, un_pw)
    if spatial_tiles:
        idxss = [find_idxs(url_func(product_name, x), selection, session) for x in tqdm.tqdm(tiles, delay = 20)]
        urls = [create_url(url_func(product_name, x), idxs, variables, request_dims = request_dims) for x, idxs in zip(tiles, idxss)]
    else:
        idxs = find_idxs(url_func(product_name, tiles[0]), selection, session)
        urls = [create_url(url_func(product_name, x), idxs, variables, request_dims = request_dims) for x in tiles]

    # Download data.
    files = download_urls(urls, "", session, fps = fps, parallel = parallel)

    # Merge spatial tiles.
    if spatial_tiles:
        dss = [process_ds(xr.open_dataset(x, decode_coords = "all"), coords, variables, crs = data_source_crs) for x in files]
        ds = rioxarray.merge.merge_datasets(dss)
    else:
        ds = process_ds(xr.open_mfdataset(files, decode_coords = "all"), coords, variables, crs = data_source_crs)

    # Reproject if necessary.
    if ds.rio.crs.to_epsg() != 4326:
        ds = ds.rio.reproject(rasterio.crs.CRS.from_epsg(4326))

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    if isinstance(timedelta, np.timedelta64):
        ds["time"] = ds["time"] + timedelta

    # Remove unnecessary coordinates.
    ds = ds.drop_vars([x for x in ds.coords if x not in ["x", "y", "time", "spatial_ref"]])

    # Save final output.
    fp = os.path.join(folder, f"{product_name}.nc")
    ds.attrs = {}
    ds = save_ds(ds, fp, decode_coords = "all")

    # Remove temporary files.
    for x in fps:
        os.remove(x)

    return ds

def find_idxs(base_url, selection, session):
    def _find_idxs(ds, k, search_range):
        all_idxs = np.where((ds[k] >= search_range[0]) & (ds[k] <= search_range[1]))[0]
        return [np.min(all_idxs), np.max(all_idxs)]
    fp = tempfile.NamedTemporaryFile(suffix=".nc").name
    url_coords = base_url + urllib.parse.quote(",".join(selection.keys()))
    fp = download_url(url_coords, fp, session, waitbar = False)
    ds = xr.open_dataset(fp, decode_coords = "all")
    idxs = {k: _find_idxs(ds, k, v) for k, v in selection.items()}
    return idxs

def create_url(base_url, idxs, variables, request_dims = True):
    if request_dims:
        dims = [f"{k}[{v[0]}:{v[1]}]" for k, v in idxs.items()]
    else:
        dims = []
    varis = [f"{k}{''.join([f'[{idxs[dim][0]}:{idxs[dim][1]}]' for dim in v[0]])}" for k, v in variables.items()]
    url = base_url + urllib.parse.quote(",".join(dims + varis))
    return url

def start_session(base_url, selection, un_pw = [None, None]):
    if un_pw == [None, None]:
        warnings.filterwarnings("ignore", "password was not set. ")
    url_coords = base_url + urllib.parse.quote(",".join(selection.keys()))
    session = setup_session(*un_pw, check_url = url_coords)
    return session

def download_xarray(url, fp, coords, variables, post_processors, 
                    data_source_crs = None, timedelta = None):

    warnings.filterwarnings("ignore", category=xr.SerializationWarning)
    online_ds = xr.open_dataset(url, decode_coords="all")
    # warnings.filterwarnings("default", category=xr.SerializationWarning)

    # Define selection.
    selection = create_selection(coords, target_crs = data_source_crs)

    # Make the selection on the remote.
    online_ds = online_ds.sel({k: slice(*v) for k, v in selection.items()})

    # Rename variables and assign crs.
    online_ds = process_ds(online_ds, coords, variables, crs = data_source_crs)

    # Download the data.
    ds = save_ds(online_ds, fp.replace(".nc", "_temp.nc"), decode_coords="all")

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    if isinstance(timedelta, np.timedelta64):
        ds["time"] = ds["time"] + timedelta

    # Save final output
    ds = save_ds(ds, fp, decode_coords="all")

    os.remove(fp.replace(".nc", "_temp.nc"))

    return ds

def create_selection(coords, ds = None, target_crs = None, 
                        source_crs = CRS.from_epsg(4326)):
    """Create a dictionary that can be given to `xr.Dataset.sel`.

    Parameters
    ----------
    coords : dict
        Dictionary describing the different dimensions over which to select. Possible keys
        are "x" for latitude, "y" for longitude and "t" for time, but other selections
        keys are also allowed (e.g. so select a band). Values are tuples with the first
        value the respective dimension names in the `ds` and the second value the selector.
    ds : xr.Dataset, optional
        Dataset on which the selection will be applied, can be supplied to check if 
        dimensions are sorted ascending or descending and to derive `target_crs`, by default None.
    target_crs : rasterio.crs.CRS, optional
        crs of the dataset on which the selection will be applied, by default None.
    source_crs : rasterio.crs.CRS, optional
        crs of the `x` and `y` limits in `coords`, by default `epsg:4326`.

    Returns
    -------
    dict
        Dimension names with slices to apply to each dimension.
    """
    selection = {}

    if not isinstance(target_crs, type(None)):
        coords["x"][1], coords["y"][1] = rasterio.warp.transform(source_crs,
                                                target_crs,
                                                coords["x"][1], coords["y"][1])
    elif not isinstance(ds, type(None)):
        if ds.rio.crs != source_crs:
            coords["x"][1], coords["y"][1] = rasterio.warp.transform(source_crs,
                                        ds.rio.crs,
                                        coords["x"][1], coords["y"][1]) 

    if "t" in coords.keys():
        coords["t"][1] = [np.datetime64(t) for t in coords["t"][1]]

    for name, lim in coords.values():
        # If `ds` is given, checks if the dimensions are sorted from high to low
        # or from low to high and adjusts the limit accordingly.
        if isinstance(ds, xr.Dataset) and hasattr(lim, "__iter__"):
            lim = {False: lim, True: lim[::-1]}[bool(ds[name][0] - ds[name][1] > 0)]
        # Create a slice object if `lim` is a list or tuple.
        # if hasattr(lim, "__iter__"):
        #     lim = slice(*lim)
        selection[name] = lim

    return selection