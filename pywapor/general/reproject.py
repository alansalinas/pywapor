import os
import rasterio
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT
import xarray as xr
import numpy as np
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
    ds_match = ds_match.assign_coords({"x": example_ds.x, "y": example_ds.y})

    label = f"Using `reproject_bulk` on {os.path.split(src_ds.encoding['source'])[-1]}:{list(src_ds.data_vars)[0]} ({spatial_interp})."
    ds_match = save_ds(ds_match, dst_path, encoding = "initiate", label = label)

    return ds_match

def reproject(src_ds, example_ds, dst_path, spatial_interp = "nearest", 
                max_bytes = 2e9, min_times = 10):

    test_ds = [src_ds, example_ds][np.argmax([src_ds.nbytes, example_ds.nbytes])]

    reproj = choose_reprojecter(test_ds, max_bytes = max_bytes, min_times = min_times)

    ds = reproj(src_ds, example_ds, dst_path, spatial_interp = spatial_interp)
    return ds

def reproject_chunk(src_ds, example_ds, dst_path, spatial_interp = "nearest"):

    das = list()
    variables = dict()
    ncs = list()

    new_src_ds = src_ds.sortby("y", ascending = False)
    new_src_ds = new_src_ds.rio.write_transform(new_src_ds.rio.transform(recalc=True))

    if not "source" in src_ds.encoding.keys() or not new_src_ds.identical(src_ds):
        src_path = src_ds.encoding.get("source", dst_path).replace(".nc", "_fixed.nc")
        src_ds = save_ds(new_src_ds, src_path, encoding = "initiate", label = "Correcting src_ds.")
        ncs.append(src_path)
    else:
        src_path = src_ds.encoding["source"]

    resampling = {'nearest': 0,
                    'bilinear': 1,
                    'cubic': 2,
                    'cubic_spline': 3,
                    'lanczos': 4,
                    'average': 5,
                    'mode': 6}[spatial_interp]

    example_ds = example_ds.sortby("y", ascending=False)

    vrt_options = {
        'resampling': resampling,
        'crs': example_ds.rio.crs,
        'transform': example_ds.rio.transform(recalc=True),
        'height': example_ds.y.size,
        'width': example_ds.x.size,
    }

    for var in src_ds.data_vars:

        part_path = dst_path.replace(".nc", f"_{var}_temp.nc")

        with rasterio.open(f'netcdf:{src_path}:{var}') as src:
            with WarpedVRT(src, **vrt_options) as vrt:
                rio_shutil.copy(vrt, part_path, driver='netcdf', creation_options = {"COMPRESS": "DEFLATE"})

        ds_part = open_ds(part_path, decode_times = False)

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
    ds = process_ds(ds, {"x": ["lon", None], "y": ["lat", None]}, variables)
    ds = ds.assign_coords({"x": example_ds.x, "y": example_ds.y})

    label = f"Using `reproject_chunk` on {os.path.split(src_ds.encoding['source'])[-1]}:{list(src_ds.data_vars)[0]} ({spatial_interp})."
    ds = save_ds(ds, dst_path, encoding = "initiate", label = label)

    for nc in ncs:
        os.remove(nc)

    return ds

if __name__ == "__main__":

    ...