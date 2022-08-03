import pywapor
from sentinelsat import SentinelAPI
import os
from datetime import datetime as dt
import tqdm
import shutil
import pywapor
from pywapor.general.processing_functions import save_ds, open_ds, create_wkt, unpack, transform_bb
from pywapor.general.logger import log
import xarray as xr
import rioxarray.merge
import tqdm
import rasterio.crs

def download(folder, latlim, lonlim, timelim, search_kwargs, node_filter = None):

    if isinstance(timelim[0], str):
        timelim[0] = dt.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = dt.strptime(timelim[1], "%Y-%m-%d")

    un, pw = pywapor.collect.accounts.get('SENTINEL')
    api = SentinelAPI(un, pw, 'https://apihub.copernicus.eu/apihub')
    
    footprint = create_wkt(latlim, lonlim)

    products = api.query(footprint, date = tuple(timelim), **search_kwargs)

    out = api.download_all(products, folder, nodefilter = node_filter)

    if isinstance(node_filter, type(None)):
        scenes = [x["path"] for x in out[0].values()]
    else:
        scenes = [os.path.join(folder, x["node_path"][2:]) for x in out[0].values()]
    
    return scenes

def process_sentinel(scenes, variables, specific_processor, time_func, final_fn, bb = None):

    example_ds = None
    dss1 = dict()

    for scene_folder in tqdm.tqdm(scenes):
        
        folder, fn = os.path.split(scene_folder)

        ext = os.path.splitext(scene_folder)[-1]
        
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

        if ext == ".zip":
            scene_folder = unpack(fn, folder)
            remove_folder = True
        else:
            scene_folder = scene_folder
            remove_folder = False

        ds = specific_processor(scene_folder, variables)

        target_crs = rasterio.crs.CRS.from_epsg(4326)

        # NOTE: see https://github.com/corteva/rioxarray/issues/545
        ds = ds.sortby("y", ascending = False)

        # Clip and pad to bounding-box
        if isinstance(example_ds, type(None)):
            if not isinstance(bb, type(None)):
                if ds.rio.crs != target_crs:
                    bb = transform_bb(target_crs, ds.rio.crs, bb)
                ds = ds.rio.clip_box(*bb)
                ds = ds.rio.pad_box(*bb)
            ds = save_ds(ds, os.path.join(scene_folder, "temp.nc")) # NOTE saving because otherwise rio.reproject bugs.
            ds = ds.rio.reproject(target_crs)
            example_ds = ds
        else:
            ds = ds.rio.reproject_match(example_ds)
            ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

        dtime = time_func(fn)
        ds = ds.expand_dims({"time": 1})
        ds = ds.assign_coords({"time":[dtime]})

        # Save to netcdf
        ds = save_ds(ds, fp)

        if dtime in dss1.keys():
            dss1[dtime].append(ds)
        else:
            dss1[dtime] = [ds]

        if remove_folder:
            shutil.rmtree(scene_folder)

    # Merge spatially.
    dss = [rioxarray.merge.merge_datasets(dss0) for dss0 in dss1.values()]

    # Merge temporally.
    ds = xr.merge(dss)

    # Define output path.
    fp = os.path.join(folder, final_fn)
    if os.path.isfile(fp):
        os.remove(fp)

    # Define encoding.
    encoding = {v: {"zlib": True, "dtype": "float32"} for v in list(ds.data_vars)}
    encoding["time"] = {"dtype": "float64"}

    # Save final netcdf.
    ds = save_ds(ds, fp, encoding = encoding)

    # Remove intermediate files.
    for dss0 in dss1.values():
        for x in dss0:
            os.remove(x.encoding["source"])

    return ds

# if __name__ == "__main__":

    # folder = r"/Users/hmcoerver/On My Mac/sentinel_dl_test/SENTINEL3"
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    # # timelim = ["2021-07-01", "2021-07-11"]
    # timelim = ["2022-06-01", "2022-06-11"]
    # product_name = 'SL_2_LST___'
    # node_filter = None

    # search_kwargs = {
    #                     "platformname": "Sentinel-3",
    #                     "producttype": product_name,
    #                     "limit": 10,
    #                     }

    # scenes = download(folder, latlim, lonlim, timelim, search_kwargs, node_filter = None)