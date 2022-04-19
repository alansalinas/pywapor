import os
from pywapor.collect.protocol import cog

def default_vars(product_name, req_vars):
    variables = {
        'GLOBCOVER': {
                "Band1": [("lat", "lon"), "lulc"],
                "crs": [(), "spatial_ref"],
                    },
    }

    req_dl_vars = {
        "GLOBCOVER": {
            "lulc": ["Band1", "crs"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}
    
    return out

def default_post_processors(product_name, req_vars):
    post_processors = {
        'GLOBCOVER': {"lulc": []},
    }

    out = [val for key, sublist in post_processors[product_name].items() for val in sublist if key in req_vars]
    if "_" in post_processors[product_name].keys():
        out += post_processors[product_name]["_"]

    return out

def url_func(product_name):
    return r"http://due.esrin.esa.int/files/GLOBCOVER_L4_200901_200912_V2.3.color.tif"

def download(folder, latlim, lonlim, product_name, req_vars = ["lulc"],
                variables = None, post_processors = None, **kwargs):
    
    folder = os.path.join(folder, "GLOBCOVER")
    
    coords = {"x": ("lon", lonlim), "y": ("lat", latlim)}

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)

    ds = cog.download(folder, product_name, coords, variables, 
                        post_processors, url_func)
    
    return ds

if __name__ == "__main__":

    product_name = 'GLOBCOVER'

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]

    ds = download(folder, latlim, lonlim, product_name,
                variables = None, post_processors = None)
