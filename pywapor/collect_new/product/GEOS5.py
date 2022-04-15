from pywapor.collect_new.protocol import opendap
from pywapor.collect_new.protocol.projections import get_crss
import os

def default_vars(product_name, req_vars = ["t_air_i", "u2m_i", "v2m_i",
                    "qv_i", "wv_i", "p_air_i", "p_air_0_i"]):

    variables = {
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

    req_dl_vars = {
        "inst3_2d_asm_Nx": {
            "t_air_i": ["t2m"],
            "u2m_i": ["u2m"],
            "v2m_i": ["v2m"],
            "qv_i": ["qv2m"],
            "wv_i": ["tqv"],
            "p_air_i": ["ps"],
            "p_air_0_i": ["slp"],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars = ["t_air_i", "u2m_i", "v2m_i",
                    "qv_i", "wv_i", "p_air_i", "p_air_0_i"]):

    post_processors = {
        "inst3_2d_asm_Nx": {
            "t_air_i": [],
            "u2m_i": [],
            "v2m_i": [],
            "qv_i": [],
            "wv_i": [],
            "p_air_i": [],
            "p_air_0_i": [],
        }
    }

    out = [val for key, sublist in post_processors[product_name].items() for val in sublist if key in req_vars]
    if "_" in post_processors[product_name].keys():
        out += post_processors[product_name]["_"]

    return out

def download(folder, latlim, lonlim, timelim, product_name, req_vars = ["t_air_i", "u2m_i", "v2m_i",
            "qv_i", "wv_i", "p_air_i", "p_air_0_i"], variables = None, post_processors = None):

    folder = os.path.join(folder, "GEOS5")

    coords = {"x": ["lon", lonlim], "y": ["lat", latlim], "t": ["time", timelim]}

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars = req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars = req_vars)

    data_source_crs = get_crss("WGS84")

    url = f"https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/{product_name}"

    fp = os.path.join(folder, f"{product_name}.nc")

    ds = opendap.download_xarray(url, fp, coords, 
                                variables, post_processors, data_source_crs = data_source_crs)
    
    return ds

if __name__ == "__main__":

    import datetime

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 11)]

    # GEOS5.
    product_name = "inst3_2d_asm_Nx"

    variables = None
    post_processors = None

    req_vars = ["t_air_i", "u2m_i", "v2m_i", "qv_i", "wv_i", "p_air_i", "p_air_0_i"]

    ds = download(folder, latlim, lonlim, timelim, product_name, req_vars = req_vars)
    print(ds.rio.crs, ds.rio.grid_mapping)