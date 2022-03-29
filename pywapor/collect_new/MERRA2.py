import datetime
import os
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
import xarray as xr
from functools import partial
import pywapor.collect.accounts as accounts
from pywapor.general.logger import log
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap
from pywapor.general.processing_functions import save_ds, create_selection

def download(folder, product_name, latlim, lonlim, timelim, variables = None,
                post_processors = None):
    
    # Define coordinate names.
    coords = {"x": "lon", "y": "lat", "t": "time"}

    # Get username and password
    un_pw = accounts.get("NASA")

    # Create array of daily dates.
    dates = pd.date_range(timelim[0], timelim[1], freq="D")

    # Define which variables to download.
    if isinstance(variables, type(None)):
        variables = default_vars(product_name)

    # Define which post processors to apply.
    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name)

    # Create selection object.
    selection = create_selection(latlim, lonlim, timelim, coords)

    # Start download session.
    idxs, session = opendap.start_session(get_url(product_name, dates[0]), selection, un_pw)

    # Prepare urls and filepaths.
    urls = [opendap.create_url(get_url(product_name, date), idxs, variables) for date in dates]
    fps = [os.path.join(folder, f"{product_name}_{date.strftime('%Y%m%d')}.nc") for date in dates]

    # Set some keyword arguments on the opendap.download_url function.
    dler = partial(opendap.download_url, session = session, waitbar = False, return_fp = True)

    # Download the urls in parallel.
    n_jobs = 8
    backend = "loky"
    files = Parallel(n_jobs=n_jobs, backend = backend)(delayed(dler)(*x) for x in tqdm(zip(urls, fps)))

    # Create the output xr.Dataset.
    ds = xr.open_mfdataset(files)
    ds = opendap.process_ds(ds, coords, variables, crs = get_crss("WGS84"))
    ds.attrs = {}

    # Apply product specific functions.
    for func in post_processors:
        ds = func(ds)

    # Save the complete netcdf.
    fp = os.path.join(folder, f"{product_name}.nc")
    ds = save_ds(ds, fp, decode_coords="all")

    return ds

def default_vars(product_name):
    vars = {
        "M2I1NXASM.5.12.4": {
                    "T2M": [("time", "lat", "lon"), "t_air_i"],
                    "U2M": [("time", "lat", "lon"), "u2m_i"],
                    "V2M": [("time", "lat", "lon"), "v2m_i"],
                    "QV2M": [("time", "lat", "lon"), "qv_i"],
                    "TQV": [("time", "lat", "lon"), "wv_i"],
                    "PS": [("time", "lat", "lon"), "p_air_i"],
                    "SLP": [("time", "lat", "lon"), "p_air_0_i"],
                        },
    }
    return vars[product_name]

def default_post_processors(product_name):
    post_processors = {
        "M2I1NXASM.5.12.4": [],
    }
    return post_processors[product_name]

def get_url(product_name, date):
    url = f"https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/{product_name}/{date.year}/{date.month:02d}/MERRA2_401.inst1_2d_asm_Nx.{date.strftime('%Y%m')}{date.day:02d}.nc4.nc4?"
    return url

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads/merra2"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 15)]

    variables = None
    post_processors = None

    product_name = "M2I1NXASM.5.12.4"

    download(folder, product_name, latlim, lonlim, timelim)