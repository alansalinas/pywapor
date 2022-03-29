import xarray as xr
import numpy as np
import os
import tqdm
import tempfile
from pydap.cas.urs import setup_session
import urllib.parse
from pywapor.general.processing_functions import save_ds, domain_overlaps_domain
import warnings

def download_url(url, fp, session, waitbar = True, return_fp = False):
    file_object = session.get(url, stream = True)
    file_object.raise_for_status()
    
    folder = os.path.split(fp)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    if os.path.isfile(fp):
        os.remove(fp)

    if waitbar:
        wb = tqdm.tqdm(unit='Bytes', unit_scale=True, position = 0)

    temp_fp = fp.replace(".nc", "_temp")

    with open(temp_fp, 'wb') as z:
        for data in file_object.iter_content(chunk_size=1024):
            size = z.write(data)
            if waitbar:
                wb.update(size)

    os.rename(temp_fp, fp)

    if return_fp:
        return fp
    else:
        ds = xr.open_dataset(fp, decode_coords = "all")
        return ds

def find_idxs(ds, k, search_range):
    all_idxs = np.where((ds[k] >= search_range[0]) & (ds[k] <= search_range[1]))[0]
    return [np.min(all_idxs), np.max(all_idxs)]

def create_url(base_url, idxs, variables, include_dims = True):
    if include_dims:
        dims = [f"{k}[{v[0]}:{v[1]}]" for k, v in idxs.items()]
    else:
        dims = []
    varis = [f"{k}{''.join([f'[{idxs[dim][0]}:{idxs[dim][1]}]' for dim in v[0]])}" for k, v in variables.items()]
    url = base_url + urllib.parse.quote(",".join(dims + varis))
    return url

def start_session(base_url, select, un_pw = [None, None]):
    if un_pw == [None, None]:
        warnings.filterwarnings("ignore", "password was not set. ")
    url_coords = base_url + urllib.parse.quote(",".join(select.keys()))
    session = setup_session(*un_pw, check_url = url_coords)
    fp = tempfile.NamedTemporaryFile(suffix=".nc").name
    ds = download_url(url_coords, fp, session, waitbar = False)
    idxs = {k: find_idxs(ds, k, v) for k, v in select.items()}
    return idxs, session

def process_ds(ds, coords, variables, crs = None):
    if isinstance(crs, type(None)):
        crs = ds.rio.crs
    ds = ds.rename({v:k for k,v in coords.items() if k in ["x", "y"]})
    ds = ds.rename({k: v[1] for k, v in variables.items()})
    if (ds.rio.grid_mapping not in list(ds.coords)) and ("spatial_ref" in [x[1] for x in variables.values()]):
        ds = ds.rio.write_grid_mapping("spatial_ref")
    ds = ds.rio.write_crs(crs)
    return ds

def opendap_to_xarray(store, fp, select, rename_keep_vars):
    
    # Open remote dataset.
    online_ds = xr.open_dataset(store, decode_coords="all")

    # Remove attributes.
    online_ds.attrs = {}

    # Check if selection is valid.
    checks = [domain_overlaps_domain(v, [online_ds[k].min(), online_ds[k].max()]) for k, v in select.items()]

    # Make the selection on the remote.
    if np.all(checks):
        online_ds = online_ds.sel({k: slice(*v) for k, v in select.items()})
    else:
        return None

    # Cleaning up.
    online_ds = online_ds[list(rename_keep_vars.keys())]
    online_ds = online_ds.rename_vars(rename_keep_vars)

    # Download the data.
    ds = save_ds(online_ds, fp, decode_coords="all")

    return ds