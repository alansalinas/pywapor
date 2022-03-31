from pywapor.collect_new import opendap
from pywapor.collect_new.projections import get_crss
import os

def default_vars(product_name):
    vars = {
        "inst3_2d_asm_Nx": {
                    "t2m": [("time", "lat", "lon"), "t_air_i"],
                    "u2m": [("time", "lat", "lon"), "u2m_i"],
                    "v2m": [("time", "lat", "lon"), "v2m_i"],
                    "qv2m": [("time", "lat", "lon"), "qv_i"],
                    "tqv": [("time", "lat", "lon"), "wv_i"],
                    "ps": [("time", "lat", "lon"), "p_air_i"],
                    "slp": [("time", "lat", "lon"), "p_air_0_i"],
                        },
    }
    return vars[product_name]

def default_post_processors(product_name):
    post_processors = {
        "inst3_2d_asm_Nx": [],
    }
    return post_processors[product_name]

def download(folder, latlim, lonlim, timelim, product_name,
                variables = None, post_processors = None):

    coords = {"x": "lon", "y": "lat", "t": "time"}

    if isinstance(variables, type(None)):
        variables = default_vars(product_name)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name)

    data_source_crs = get_crss("WGS84")

    url = f"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/{product_name}"

    fp = os.path.join(folder, f"{product_name}.nc")

    ds = opendap.download_xarray(url, fp, latlim, lonlim, timelim, coords, 
                                variables, post_processors, data_source_crs = data_source_crs)
    
    return ds

if __name__ == "__main__":

    import datetime
    folder = r"/Users/hmcoerver/Downloads/merra2"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 8, 15)]

    # GEOS5.
    product_name = "inst3_2d_asm_Nx"

    variables = None
    post_processors = None

    ds = download(folder, latlim, lonlim, timelim, product_name)
