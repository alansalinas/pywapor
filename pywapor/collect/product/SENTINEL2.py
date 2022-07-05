import glob
import os
import shutil
import xarray as xr
import numpy as np
import pywapor
from pywapor.collect.product.Landsat.C2L2SP import calc_albedo, calc_ndvi, mask_invalid, open_as_xr, transform_bb
from datetime import datetime as dt
import rasterio.crs
from pywapor.general.logger import log
from pywapor.general.processing_functions import save_ds, open_ds
from sentinelsat import SentinelAPI
import numpy as np
from pywapor.enhancers.apply_enhancers import apply_enhancer
import rioxarray.merge

def processor(scenes, variables, bb = None):

    example_ds = None
    dss1 = dict()

    for scene_folder in scenes:

        folder, fn = os.path.split(scene_folder)

        fp = os.path.join(folder, os.path.splitext(fn)[0] + ".nc")
        if os.path.isfile(fp):
            ds = open_ds(fp)
            dtime = ds.time.values[0]
            if dtime in dss1.keys():
                dss1[dtime].append(ds)
            else:
                dss1[dtime] = [ds]
            continue

        log.info(f"--> Processing `{fn}`.")

        if os.path.splitext(fn)[-1] == ".zip":
            scene_folder = unpack(fn, folder)

        fps = {v[1]: glob.glob(os.path.join(scene_folder, "**", "*" + k), recursive = True)[0] for k, v in variables.items()}
        ds = xr.concat([open_as_xr(fp, name) for name, fp in fps.items()], "band")

        qa_da = ds.band_data.sel(band = "qa")
        data = ds.where(ds.band != "qa", drop = True)

        valid_range = {
                "blue":     (1, 65534), # 0 = NODATA, 65535 = SATURATED
                "green":    (1, 65534),
                "red":      (1, 65534),
                "nir":      (1, 65534),
                }

        masked_data = mask_invalid(data, valid_range)

        # 0 SC_NODATA # 1 SC_SATURATED_DEFECTIVE # 2 SC_DARK_FEATURE_SHADOW
        # 3 SC_CLOUD_SHADOW # 4 SC_VEGETATION # 5 SC_NOT_VEGETATED
        # 6 SC_WATER # 7 SC_UNCLASSIFIED # 8 SC_CLOUD_MEDIUM_PROBA
        # 9 SC_CLOUD_HIGH_PROBA # 10 SC_THIN_CIRRUS # 11 SC_SNOW_ICE
        pixel_qa_flags = [0, 1, 2, 3, 7, 8, 9, 10, 11]
        keep = np.invert(qa_da.isin(pixel_qa_flags))
        masked_data = masked_data.where(keep)

        scale = 1./10000. # BOA_QUANTIFICATION_VALUE
        offset = -1000 # BOA_ADD_OFFSET
        scaled_data = (masked_data + offset) * scale

        scaled_data = scaled_data.where(scaled_data <= 1.00)
        scaled_data = scaled_data.where(scaled_data >= 0.00)

        das = list()

        if np.all([x in scaled_data.band.values for x in ["red", "nir"]]):
            das.append(calc_ndvi(scaled_data))

        weights = {
            "blue": 0.074,
            "green": 0.083,
            "red": 0.334,
            "nir": 0.356,
            "offset": 0.033,
        }

        if np.all([x in scaled_data.band.values for x in ["blue", "green", "red", "nir"]]):
            das.append(calc_albedo(scaled_data, weights = weights))

        ds = xr.merge(das)

        # Add time dimension to data arrays.
        ds = ds.expand_dims({"time": 1})
        dtime = np.datetime64(dt.strptime(fn.split("_")[2], "%Y%m%dT%H%M%S"))
        ds = ds.assign_coords({"time":[dtime]})

        target_crs = rasterio.crs.CRS.from_epsg(4326)

        # Clip and pad to bounding-box
        if isinstance(example_ds, type(None)):
            if not isinstance(bb, type(None)):
                local_bb = transform_bb(target_crs, ds.rio.crs, bb)
                ds = ds.rio.clip_box(*local_bb)
                ds = ds.rio.pad_box(*local_bb)
            ds = save_ds(ds, os.path.join(scene_folder, os.path.splitext(fn)[0], "temp.nc")) # NOTE saving because otherwise rio.reproject bugs.
            ds = ds.rio.reproject(target_crs)
            example_ds = ds
        else:
            ds = ds.rio.reproject_match(example_ds)
            ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

        # Save to netcdf
        ds = save_ds(ds, fp)

        if dtime in dss1.keys():
            dss1[dtime].append(ds)
        else:
            dss1[dtime] = [ds]

    # Merge spatially.
    dss = [rioxarray.merge.merge_datasets(dss0) for dss0 in dss1.values()]

    # Merge temporally.
    ds = xr.merge(dss, combine_attrs = "drop")

    fp = os.path.join(folder, "SENTINEL2.nc")

    if os.path.isfile(fp):
        os.remove(fp)

    encoding = {v: {"zlib": True, "dtype": "float32"} for v in list(ds.data_vars)}
    encoding["time"] = {"dtype": "float64"}
    ds = save_ds(ds, fp, encoding = encoding)

    for dss0 in dss1.values():
        for x in dss0:
            os.remove(x.encoding["source"])

    return ds

def unpack(file, folder):
    fn = os.path.splitext(file)[0]
    shutil.unpack_archive(os.path.join(folder, file), folder)
    folder = [x for x in glob.glob(os.path.join(folder, fn + "*")) if os.path.isdir(x)][0]
    return folder

def create_wkt(latlim, lonlim):
    left = lonlim[0]
    bottom = latlim[0]
    right = lonlim[1]
    top = latlim[1]
    x = f"{left} {bottom},{right} {bottom},{right} {top},{right} {bottom},{left} {bottom}"
    return "GEOMETRYCOLLECTION(POLYGON((" + x + ")))"

def default_vars(product_name, req_vars):

    variables = {
        "S2MSI2A": {
                    "_B02_20m.jp2": [(), "blue"],
                    "_B03_20m.jp2": [(), "green"],
                    "_B04_20m.jp2": [(), "red"],
                    "_B8A_20m.jp2": [(), "nir"],
                    "_SCL_20m.jp2": [(), "qa"],
                },
    }

    req_dl_vars = {
        "S2MSI2A": {
            "ndvi": ["_B04_20m.jp2", "_B8A_20m.jp2", "_SCL_20m.jp2"],
            "r0": ["_B02_20m.jp2", "_B03_20m.jp2", "_B04_20m.jp2", "_B8A_20m.jp2", "_SCL_20m.jp2"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out

def default_post_processors(product_name, req_vars):
    return {}

def download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = None, post_processors = None):

    product_folder = os.path.join(folder, "SENTINEL2")

    fn = os.path.join(product_folder, f"{product_name}.nc")

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    if isinstance(timelim[0], str):
        timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    un, pw = pywapor.collect.accounts.get('SENTINEL')
    api = SentinelAPI(un, pw, 'https://apihub.copernicus.eu/apihub')
    footprint = create_wkt(latlim, lonlim)

    products = api.query(
                        footprint,
                        date = tuple(timelim),
                        platformname = "Sentinel-2",
                        producttype = product_name,
                        cloudcoverpercentage = (0, 30),
                        limit = 10
                        )

    def node_filter(node_info):
        fn = os.path.split(node_info["node_path"])[-1]
        to_dl = list(variables.keys())
        return np.any([x in fn for x in to_dl])

    out = api.download_all(products, product_folder, nodefilter = node_filter)

    scenes = [os.path.join(product_folder, x["node_path"][2:]) for x in out[0].values()]

    bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]
    ds = processor(scenes, variables, bb = bb)

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/On My Mac/sentinel_dl_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    # timelim = ["2021-07-01", "2021-07-11"]
    timelim = ["2022-06-01", "2022-06-11"]
    product_name = 'S2MSI2A'
    # product_name = 'SL_2_LST___'

    req_vars = ["ndvi", "r0"]
    # req_vars = ["lst"]
    post_processors = None
    variables = None

    check = download(folder, latlim, lonlim, timelim, product_name, 
                req_vars, variables = None,  post_processors = None)
    
    # left = lonlim[0]
    # bottom = latlim[0]
    # right = lonlim[1]
    # top = latlim[1]
    # final_bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

    # resolution = 0.01

    # un, pw = pywapor.collect.accounts.get('SENTINEL')

    # api = SentinelAPI(un, pw, 'https://apihub.copernicus.eu/apihub')

    # footprint = create_wkt(latlim, lonlim)

    # if isinstance(timelim[0], str):
    #     timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
    #     timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    # products = api.query(footprint,
    #                      date=(timelim[0], timelim[1]),
    #                     #  platformname = 'Sentinel-3',
    #                      producttype = 'SL_2_LST___',
    #                     #  cloudcoverpercentage = (0, 30)
    #                      )



# import xarray as xr
# import glob

# nc_files = glob.glob(r"/Users/hmcoerver/On My Mac/ndvi_r0_test/S3A_SL_2_LST____20220610T083131_20220610T083431_20220611T173804_0179_086_178_2340_PS1_O_NT_004.SEN3/*_in.nc")

# ds = xr.open_mfdataset(nc_files)[["longitude_in", "latitude_in", "LST"]]
# ds = ds.set_coords(("longitude_in", "latitude_in"))
# ds = ds.rename_vars({"longitude_in": "x", "latitude_in": "y"})
# ds = ds.rename_dims({"rows": "ny", "columns": "nx"})





# test_ds = regrid_curvilinear(ds, resolution, bb = bb)




