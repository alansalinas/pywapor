import os
import rasterio
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT
import xarray as xr
import os
from pywapor.general.logger import log
from pywapor.general.processing_functions import process_ds
from pywapor.general.processing_functions import open_ds, save_ds

def choose_reprojecter(src_ds, max_bytes = 2e9, min_times = 10):

    if "time" in src_ds.dims:
        tsize = src_ds.dims["time"]
    else:
        tsize = 1

    if src_ds.nbytes < max_bytes and tsize > min_times:
        reproject = reproject_bulk
    else:
        reproject = reproject_chunk

    return reproject

def reproject_bulk(src_ds, example_ds, dst_path, spatial_interp = "nearest"):

    resampling = {'nearest': 0,
                    'bilinear': 1,
                    'cubic': 2,
                    'cubic_spline': 3,
                    'lanczos': 4,
                    'average': 5,
                    'mode': 6}[spatial_interp]
    
    ds_match = src_ds.rio.reproject_match(example_ds, resampling = resampling)
    ds_match = ds_match.assign_coords({
                                        "x": example_ds.x,
                                        "y": example_ds.y,
                                    })

    ds_match = save_ds(ds_match, dst_path, decode_coords="all")

    return ds_match

def reproject(src_ds, example_ds, dst_path, spatial_interp = "nearest", 
                max_bytes = 2e9, min_times = 10):
    reproj = choose_reprojecter(src_ds, max_bytes = max_bytes, min_times = min_times)
    log.info(f"--> Using `{reproj.__name__}` on {os.path.split(src_ds.encoding['source'])[-1]}:{list(src_ds.data_vars)[0]} ({spatial_interp}).")
    ds = reproj(src_ds, example_ds, dst_path, spatial_interp = spatial_interp)
    return ds

def reproject_chunk(src_ds, example_ds, dst_path, spatial_interp = "nearest"):

    src_path = src_ds.encoding["source"]

    # NOTE: not using encoding because it messes up the rioxarray projections
    # TODO: File bug-report.
    # encoding = {var: {"dtype": str(src_ds[var].dtype)} for var in src_ds.data_vars}

    resampling = {'nearest': 0,
                    'bilinear': 1,
                    'cubic': 2,
                    'cubic_spline': 3,
                    'lanczos': 4,
                    'average': 5,
                    'mode': 6}[spatial_interp]

    vrt_options = {
        'resampling': resampling,
        'crs': example_ds.rio.crs,
        'transform': example_ds.rio.transform(),
        'height': example_ds.y.size,
        'width': example_ds.x.size,
        'src_crs': src_ds.rio.crs,
        'src_transform': src_ds.rio.transform(),
    }

    das = list()
    variables = dict()
    ncs = list()

    for var in src_ds.data_vars:

        part_path = dst_path.replace(".nc", f"_{var}_temp.nc")

        with rasterio.open(f'netcdf:{src_path}:{var}') as src:
            with WarpedVRT(src, **vrt_options) as vrt:
                rio_shutil.copy(vrt, part_path, driver='netcdf')

        ds_part = xr.open_dataset(part_path, chunks = "auto", decode_coords="all", decode_times = False)

        if "time" in src_ds.coords:
            da = ds_part.to_array("time", name = var).assign_coords({"time": src_ds.time})
        else:
            da = ds_part[var]

        das.append(da)
        ncs.append(part_path)
        variables[var] = (None, var)

    ds = xr.merge(das)

    if "spatial_ref" not in ds.coords:
        ds = ds.assign_coords({"spatial_ref": ds.crs})

    ds = ds.drop_vars("crs")

    coords = {"x": ["lon", None], "y": ["lat", None]}

    ds = process_ds(ds, coords, variables)

    ds = ds.assign_coords({
                            "x": example_ds.x,
                            "y": example_ds.y,
                        })

    ds = save_ds(ds, dst_path, decode_coords = "all")#, encoding = encoding)

    for nc in ncs:
        os.remove(nc)

    return ds

if __name__ == "__main__":

    src_ds = xr.open_dataset(r"/Users/hmcoerver/On My Mac/create_table/MODIS/MCD43A3.061.nc", decode_coords="all")
    example_ds = xr.open_dataset(r"/Users/hmcoerver/On My Mac/create_table/MODIS/MOD13Q1.061.nc", decode_coords="all")
    dst_path = r"/Users/hmcoerver/On My Mac/create_table/test.nc"
    out = reproject_chunk(src_ds, example_ds, dst_path, spatial_interp = "nearest")

    # import tracemalloc
    # import datetime
    # import glob

    # example_path = r"/Users/hmcoerver/On My Mac/create_table/MODIS/MOD13Q1.061.nc"
    # example_ds = xr.open_dataset(example_path, decode_coords="all")

    # folder = r"/Users/hmcoerver/On My Mac/create_table/MERRA2"
    # src_paths = glob.glob(os.path.join(folder, "*.nc"))

    # for src_path in src_paths:

    #     _, fn = os.path.split(src_path)
    #     dst_path = os.path.join(folder, fn)

    #     reproject = choose_reprojecter(src_path, max_bytes = 2e9, min_times = 10)

    #     print(f"--> Using `{reproject.__name__}` for {fn}.")

    #     t1 = datetime.datetime.now()

    #     dst_path = dst_path.replace(".nc", f"_{reproject.__name__.split('_')[-1]}.nc")
    #     if os.path.isfile(dst_path):
    #         dst_path = dst_path.replace(".nc", "_.nc")

    #     tracemalloc.start()
    #     ds1 = reproject(src_path, example_ds, dst_path)
    #     mem_test1 = tracemalloc.get_traced_memory()
    #     print(f"--> Memory usage: {mem_test1[1]-mem_test1[0]}")
    #     tracemalloc.stop()

    #     t2 = datetime.datetime.now()
    #     print(f"--> Time: {t2-t1}\n")
