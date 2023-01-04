import pywapor
import inspect
import datetime
import os
import warnings
import rasterio
import glob
import numpy as np
from pywapor.general.logger import adjust_logger, log

def get_source_details(mod):
    # Parse the `default_vars`-function code content.
    lines = inspect.getsourcelines(mod.default_vars)

    # Read the lines from `default_vars` to substract the availabel products and variables.
    idx1 = [i for i, x in enumerate(lines[0]) if "variables = " in x][0]
    idx2 = [i - 1 for i, x in enumerate(lines[0]) if "req_dl_vars = " in x][0]
    ldict = {}
    exec("".join(lines[0][idx1:idx2]).replace("  ", ""), globals(), ldict)
    variables = ldict["variables"]

    # Read the lines from `default_vars` to subtract the `req_dl_vars`
    idx1 = [i for i, x in enumerate(lines[0]) if "req_dl_vars = " in x][0]
    idx2 = [i - 1 for i, x in enumerate(lines[0]) if "out = " in x][0]
    ldict = {}
    exec("".join(lines[0][idx1:idx2]).replace("  ", ""), globals(), ldict)
    req_dl_vars = ldict["req_dl_vars"]

    # List the available products for this `source`.
    products = variables.keys()

    return req_dl_vars, products

def download_products(mod, products, req_dl_vars, args):
    dss = {}

    for product_name in products:

        req_vars = list(req_dl_vars[product_name].keys())
        
        args.update({
            "product_name": product_name,
            "req_vars": req_vars,
        })

        dss[product_name] = mod.download(**args)

    return dss

def test_download(product_name, mod, workdir, timelim, latlim = [29.4, 29.7], lonlim = [30.7, 31.0], empty_folder = True, create_subfolder = True):
    
    if create_subfolder:
        folder = os.path.join(workdir, product_name)
    else:
        folder = workdir

    if not os.path.isdir(folder):
        os.makedirs(folder)
    if empty_folder:
        for f in glob.glob(os.path.join(folder, "*")):
            os.remove(f)

    adjust_logger(True, folder, "INFO")
    log.info(f"## {product_name}")

    req_dl_vars, products = get_source_details(mod)

    args = {
        "folder": folder,
        "latlim": latlim,
        "lonlim": lonlim,
        "timelim": timelim, 
    }
    dss = download_products(mod, products, req_dl_vars, args)

    return dss

def has_geotransform(ds):
    varis = ds.data_vars
    for var in varis:
        with warnings.catch_warnings(record=True) as w:
            _ = rasterio.open(f"netcdf:{ds.encoding['source']}:{var}")
            if len(w) > 0:
                for warning in w:
                    no_geot = "Dataset has no geotransform, gcps, or rpcs." in str(warning.message)
                    if no_geot:
                        return False
    return True

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

if __name__ == "__main__":

    # workdir = r"/Users/hmcoerver/Local/20220325_20220415_test_data"
    workdir = r"/Users/hmcoerver/Local/collect_tests"

    # overwrite_timelim = [datetime.date(2022, 3, 25), datetime.date(2022, 4, 15)]
    overwrite_timelim = None

    sources = {
        'GEOS5':        [pywapor.collect.product.GEOS5,     [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]], # opendap.download_xarray
        'STATICS':      [pywapor.collect.product.STATICS,   None], # cog.download
        'MODIS':        [pywapor.collect.product.MODIS,     [datetime.date(2019, 3, 1), datetime.date(2019, 4, 1)]], # opendap.download
        'MERRA2':       [pywapor.collect.product.MERRA2,    [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]], # opendap.download
        'GLOBCOVER':    [pywapor.collect.product.GLOBCOVER, None], # cog.download
        'CHIRPS':       [pywapor.collect.product.CHIRPS,    [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]], # opendap.download
        'SRTM':         [pywapor.collect.product.SRTM,      None], # opendap.download
        'PROBAV':       [pywapor.collect.product.PROBAV,    [datetime.date(2021, 7, 1), datetime.date(2021, 7, 11)]],
        'ERA5':         [pywapor.collect.product.ERA5,      [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]], # cds.download
        'SENTINEL2':    [pywapor.collect.product.SENTINEL2, [datetime.date(2022, 3, 1), datetime.date(2022, 3, 9)]],    
        'SENTINEL3':    [pywapor.collect.product.SENTINEL3, [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]],
        'VIIRSL1':      [pywapor.collect.product.VIIRSL1,   [datetime.date(2022, 3, 1), datetime.date(2022, 3, 3)]],
        'COPERNICUS':   [pywapor.collect.product.COPERNICUS, None], # cog.download
    }

    for product_name, (mod, timelim) in sources.items():
        if not isinstance(overwrite_timelim, type(None)):
            timelim = overwrite_timelim
        latlim = [29.4, 29.7]
        lonlim = [30.7, 31.0]
        dss = test_download(product_name, mod, workdir, timelim, latlim = latlim, lonlim = lonlim, empty_folder = False, create_subfolder = False)
        for ds in dss.values():
            assert ds.rio.crs.to_epsg() == 4326
            assert "spatial_ref" in ds.coords
            assert strictly_increasing(ds.x.values)
            assert strictly_decreasing(ds.y.values)
            if "time" in ds.dims:
                assert strictly_increasing(ds.time.values)
            assert has_geotransform(ds)
            assert np.all([int(ds[var].notnull().sum().values) > 0 for var in ds.data_vars])
            # for var in ds.data_vars:
            #     print(product_name, ds[var].encoding)
            #     break
