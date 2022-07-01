import glob
import os
import shutil
import xarray as xr
import numpy as np
from pywapor.collect.product.Landsat.C2L2SP import calc_albedo, calc_ndvi, mask_invalid, open_as_xr, transform_bb
from datetime import datetime as dt
import rasterio
from pywapor.general.logger import log
from pywapor.general.processing_functions import save_ds
from sentinelsat import SentinelAPI

def main(folder, final_bb = None):
    # Look for tar-files in the input folder.
    all_files = glob.glob(os.path.join(folder, "*.zip"))

    # Filter zip-files.
    files = [os.path.split(file)[-1] for file in all_files if "MSIL2A" in file]

    dss = list()

    example_ds = None

    # Loop over the zip-files.
    for file in files:
        
        log.info(f"--> Processing `{file}`.")

        # Unpack tar-file.
        scene_folder = unpack(file, folder)

        to_open = ["blue", "red", "green", "nir"]

        resolution = "60m"
        
        search_paths = {
            "blue": f"GRANULE/L2A_*/IMG_DATA/R{resolution}/*_B02*_{resolution}.jp2",
            "green": f"GRANULE/L2A_*/IMG_DATA/R{resolution}/*_B03*_{resolution}.jp2",
            "red": f"GRANULE/L2A_*/IMG_DATA/R{resolution}/*_B04*_{resolution}.jp2",
            "nir": f"GRANULE/L2A_*/IMG_DATA/R{resolution}/*_B*8*_{resolution}.jp2",
        }

        fps = {k: glob.glob(os.path.join(scene_folder, v))[0] for k, v in search_paths.items() if k in to_open}
        
        data = xr.concat([open_as_xr(fp, name) for name, fp in fps.items()], "band")

        valid_range = {
                "blue":     (1, 65534), # 0 = NODATA, 65535 = SATURATED
                "green":    (1, 65534),
                "red":      (1, 65534),
                "nir":      (1, 65534),
                }

        masked_data = mask_invalid(data, valid_range)

        qa_search_path = f"GRANULE/L2A_*/IMG_DATA/R{resolution}/*_SCL_{resolution}.jp2"
        qa_fp = glob.glob(os.path.join(scene_folder, qa_search_path))

        if len(qa_fp) == 1:

            qa_data = xr.concat([open_as_xr(fp, name) for name, fp in {"pixel_qa": qa_fp[0]}.items()], "band")
        
            pixel_qa_flags = [0, 1, 2, 3, 7, 8, 9, 10, 11]

            keep = np.invert(qa_data.band_data.isin(pixel_qa_flags))

            masked_data = masked_data.where(keep.isel(band = 0))

        scale = 1./10000. # BOA_QUANTIFICATION_VALUE
        offset = -1000 # BOA_ADD_OFFSET

        scaled_data = (masked_data + offset) * scale

        scaled_data = scaled_data.where(scaled_data <= 1.00)
        scaled_data = scaled_data.where(scaled_data >= 0.00)

        ndvi = calc_ndvi(scaled_data)

        weights = {
            "blue": 0.074,
            "green": 0.083,
            "red": 0.334,
            "nir": 0.356,
            "offset": 0.033,
        }

        albedo = calc_albedo(scaled_data, weights = weights)

        # Merge the three variables into one dataset.
        ds = xr.merge([ndvi, albedo])

        # Add time dimension to data arrays.
        ds = ds.expand_dims({"time": 1})

        dtime = np.datetime64(dt.strptime(file.split("_")[2], "%Y%m%dT%H%M%S"))
        ds = ds.assign_coords({"time":[dtime]})

        target_crs = rasterio.crs.CRS.from_epsg(4326)

        # Clip and pad to bounding-box
        if isinstance(example_ds, type(None)):
            if not isinstance(final_bb, type(None)):
                bb = transform_bb(target_crs, ds.rio.crs, final_bb)
                ds = ds.rio.clip_box(*bb)
                ds = ds.rio.pad_box(*bb)
            ds = save_ds(ds, os.path.join(scene_folder, os.path.splitext(file)[0], "temp.nc")) # NOTE saving because otherwise rio.reproject bugs.
            ds = ds.rio.reproject(target_crs)
            example_ds = ds
        else:
            ds = ds.rio.reproject_match(example_ds)
            ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

        # Save to netcdf
        fp = os.path.join(folder, os.path.splitext(file)[0] + ".nc")
        ds = save_ds(ds, fp)

        dss.append(ds)

        shutil.rmtree(scene_folder)

    ds = xr.merge(dss, combine_attrs = "drop")

    fp = os.path.join(folder, "SENTINEL2.nc")

    if os.path.isfile(fp):
        os.remove(fp)

    encoding = {v: {"zlib": True, "dtype": "float32"} for v in list(ds.data_vars)}
    encoding["time"] = {"dtype": "float64"}

    ds = save_ds(ds, fp, encoding = encoding)

    for x in dss:
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

def download(folder, latlim, lonlim, timelim, product_name, req_vars = None, post_processors = None):

    un, pw = ("","")

    api = SentinelAPI(un, pw, 'https://apihub.copernicus.eu/apihub')

    footprint = create_wkt(latlim, lonlim)

    if isinstance(timelim[0], str):
        timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    if product_name != 'S2MSI2A':
        log.warning(f"--> {product_name} is not supported, using 'S2MSI2A' instead.")
        product_name = 'S2MSI2A'

    products = api.query(footprint,
                         date=(timelim[0], timelim[1]),
                         platformname = 'Sentinel-2',
                         producttype = 'S2MSI2A',
                         cloudcoverpercentage = (0, 30))

    def node_filter(node_info):
        fn = os.path.split(node_info["node_path"])[-1]
        to_dl = ["_B02_20m.jp2", "_B03_20m.jp2", "_B04_20m.jp2", 
                    "_B8A_20m.jp2", "_SCL_20m.jp2"]
        return np.any([x in fn for x in to_dl])

    out = api.download_all(products, folder, nodefilter=node_filter)

    return out

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/On My Mac/ndvi_r0_test"

    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2021-07-01", "2021-07-11"]
    product_name = 'S2MSI2A'

    req_vars = None
    post_processors = None
    
    left = lonlim[0]
    bottom = latlim[0]
    right = lonlim[1]
    top = latlim[1]
    final_bb = [left, bottom, right, top]

# 0 SC_NODATA
# 1 SC_SATURATED_DEFECTIVE
# 2 SC_DARK_FEATURE_SHADOW
# 3 SC_CLOUD_SHADOW
# 4 SC_VEGETATION
# 5 SC_NOT_VEGETATED
# 6 SC_WATER
# 7 SC_UNCLASSIFIED
# 8 SC_CLOUD_MEDIUM_PROBA
# 9 SC_CLOUD_HIGH_PROBA
# 10 SC_THIN_CIRRUS
# 11 SC_SNOW_ICE