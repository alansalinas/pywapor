import datetime
import os
from pywapor.general.logger import log
from shapely.geometry import shape
from pywapor.collect_new.projections import get_crss
import pywapor.collect_new.opendap as opendap
from pywapor.general.processing_functions import save_ds, create_selection, open_ds

def download(folder, latlim, lonlim, timelim, variables = None):

    coords = {"x": "longitude", "y": "latitude", "t": "time"}

    # Define output name.
    fn = f"CHIRPS.nc"
    fp = os.path.join(folder, fn)
    fp_temp = fp.replace(".nc", "_temp.nc")

    log.info(f"--> Downloading {fn}.")

    # Create selection object.
    selection = create_selection(latlim, lonlim, timelim, coords)

    # Define which variables to download.
    if isinstance(variables, type(None)):
        variables = default_vars()

    if os.path.isfile(fp):
        # Open existing dataset.
        ds = open_ds(fp, decode_coords = "all")
    else:
        # Define CHIRPS tile URL.
        base_url = get_url()

        # Start OPeNDAP session.
        idxs, session = opendap.start_session(base_url, selection)
        
        # Make data request URL.
        url_data = opendap.create_url(base_url, idxs, variables, include_dims = False)
        
        # Download data.
        ds_data = opendap.download_url(url_data, fp_temp, session, coords)

        # Rename (spatial) variables.
        ds_data = opendap.process_ds(ds_data, coords, variables, crs = get_crss("WGS84"))

        ds = save_ds(ds_data, fp, decode_coords = "all")

        os.remove(fp_temp)

    return ds

def default_vars():
    variables =  {
        "precip": [("time", "latitude", "longitude"), "p"],
            }
    return variables

def get_url():
    url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/chirps20GlobalDailyP05.nc?"
    return url

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Downloads"

    product_name = "CHIRPS"
    latlim = [36.9, 43.7]
    lonlim = [5.2, 17.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 1)]

    ds = download(folder, latlim, lonlim, timelim)