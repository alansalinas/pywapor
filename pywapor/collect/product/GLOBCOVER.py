import os
from pywapor.collect.protocol import cog
from pywapor.general.processing_functions import open_ds
from functools import partial
from pywapor.enhancers import lulc

def default_vars(product_name, req_vars):
    variables = {
        '2009_V2.3_Global': {
                "Band1": [("lat", "lon"), "lulc"],
                "crs": [(), "spatial_ref"],
                    },
    }

    req_dl_vars = {
        "2009_V2.3_Global": {
            "lulc": ["Band1", "crs"],
            "rs_min": ["Band1", "crs"],
            "z_obst_max": ["Band1", "crs"], 
            "land_mask": ["Band1", "crs"],
            "lue_max": ["Band1", "crs"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def remove_var(ds, var):
    return ds.drop_vars([var])

def default_post_processors(product_name, req_vars):

    post_processors = {
        '2009_V2.3_Global': {
            "lulc": [],
            "rs_min": [partial(lulc.lulc_to_x, in_var = "lulc", out_var = "rs_min", 
                        convertor = lulc.globcover_to_rs_min())],
            "z_obst_max": [partial(lulc.lulc_to_x, in_var = "lulc", out_var = "z_obst_max", 
                        convertor = lulc.globcover_to_z_obst_max())],
            "land_mask": [partial(lulc.lulc_to_x, in_var = "lulc", out_var = "land_mask", 
                        convertor = lulc.globcover_to_land_mask())],
            "lue_max": [partial(lulc.lulc_to_x, in_var = "lulc", out_var = "lue_max", 
                        convertor = lulc.globcover_to_lue_max())],
            },
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def url_func(product_name):
    return r"http://due.esrin.esa.int/files/GLOBCOVER_L4_200901_200912_V2.3.color.tif"

def download(folder, latlim, lonlim, product_name, req_vars = ["lulc"],
                variables = None, post_processors = None, **kwargs):
    
    folder = os.path.join(folder, "GLOBCOVER")

    fn = os.path.join(folder, f"{product_name}.nc")
    if os.path.isfile(fn):
        return open_ds(fn, "all")

    coords = {"x": ("lon", lonlim), "y": ("lat", latlim)}

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    ds = cog.download(folder, product_name, coords, variables, 
                        post_processors, url_func)
    
    return ds

# if __name__ == "__main__":

#     product_name = '2009_V2.3_Global'

#     folder = r"/Users/hmcoerver/Downloads/pywapor_test"
#     # latlim = [26.9, 33.7]
#     # lonlim = [25.2, 37.2]
#     latlim = [28.9, 29.7]
#     lonlim = [30.2, 31.2]

#     ds = download(folder, latlim, lonlim, product_name,
#                 variables = None, post_processors = None)
