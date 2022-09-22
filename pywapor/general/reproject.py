import os
import rasterio
from rasterio import shutil as rio_shutil
from rasterio.vrt import WarpedVRT
import xarray as xr
import numpy as np
import os
from pywapor.general.logger import log
from pywapor.general.performance import performance_check
from pywapor.general.processing_functions import process_ds
from pywapor.general.processing_functions import open_ds, save_ds, remove_ds, log_example_ds
import rasterio
import warnings

def get_pixel_sizes(dss):
    # Check CRSs of datasets.
    crss = [v.rio.crs.to_epsg() for v in dss]
    # Count occurence of the different CRSs.
    uniqs, counts = np.unique(crss, return_counts=True)
    # Pick the dominant CRS.
    crs = uniqs[np.argmax(counts)]
    # Reproject to common CRS.
    dss2 = [(ds, ds.rio.reproject(rasterio.CRS.from_epsg(crs))) for ds in dss]
    return {np.abs(np.prod(ds2.rio.resolution())): ds for ds, ds2 in dss2}

def has_changed(ds):
    if "source" in ds.encoding.keys():
        changed = not ds.equals(open_ds(ds.encoding["source"]))
    else:
        changed = True
    return changed

def has_geotransform(ds):
    varis = ds.data_vars
    if not "source" in ds.encoding.keys():
        return False
    for var in varis:
        with warnings.catch_warnings(record=True) as w:
            _ = rasterio.open(f"netcdf:{ds.encoding['source']}:{var}")
            if len(w) > 0:
                for warning in w:
                    no_geot = "Dataset has no geotransform, gcps, or rpcs." in str(warning.message)
                    if no_geot:
                        return False
    return True

def align_pixels(dss, folder, spatial_interp = "nearest", example_ds = None, stack_dim = "time", fn_append = ""):

    temp_files = list()

    if isinstance(spatial_interp, str):
        spatial_interp = [spatial_interp] * len(dss)
    assert len(dss) == len(spatial_interp)

    if len(dss) == 1 and isinstance(example_ds, type(None)):
        dss1 = [dss[0]]
    else:
        if isinstance(example_ds, type(None)):
            example_ds = min(get_pixel_sizes(dss).items(), key=lambda x: x[0])[1]
            log_example_ds(example_ds)
        dss1 = list()
        for i, (spat_interp, ds_part) in enumerate(zip(spatial_interp, dss)):
            if needs_reprojecting(ds_part, example_ds):
                var_str = "_".join(ds_part.data_vars)
                dst_path = os.path.join(folder, f"{var_str}_x{i}{fn_append}.nc")
                ds_part = reproject(ds_part, example_ds, dst_path, spatial_interp = spat_interp, stack_dim = stack_dim)
                temp_files.append(ds_part)
            else:
                ds_part = ds_part.assign_coords({"y": example_ds.y, 
                                                "x": example_ds.x, 
                                                "spatial_ref": example_ds.spatial_ref})
            dss1.append(ds_part)

    return dss1, temp_files

def needs_reprojecting(ds, example_ds):
    check_crs = ds.rio.crs == example_ds.rio.crs
    check_ysize = ds.y.size == example_ds.y.size
    check_xsize = ds.x.size == example_ds.x.size
    check_y = bool((ds.y == example_ds.y).all().values)
    check_x = bool((ds.x == example_ds.x).all().values)
    return not np.all([check_crs, check_ysize, check_xsize, check_x, check_y])

def choose_reprojecter(src_ds, max_bytes = 2e9, min_times = 10, stack_dim = "time"):

    if stack_dim in src_ds.dims:
        tsize = src_ds.dims[stack_dim]
    else:
        tsize = 1

    if src_ds.nbytes < max_bytes and tsize > min_times:
        reproject = reproject_bulk
    else:
        reproject = reproject_chunk

    return reproject

def reproject_bulk(src_ds, example_ds, dst_path, spatial_interp = "nearest", **kwargs):

    resampling = {'nearest': 0,
                    'bilinear': 1,
                    'cubic': 2,
                    'cubic_spline': 3,
                    'lanczos': 4,
                    'average': 5,
                    'mode': 6}[spatial_interp]
    
    ds_match = src_ds.rio.reproject_match(example_ds, resampling = resampling)
    ds_match = ds_match.assign_coords({"x": example_ds.x, "y": example_ds.y})

    label = f"Applying `reproject_bulk` to {os.path.split(src_ds.encoding['source'])[-1]}:{list(src_ds.data_vars)[0]} ({spatial_interp})."
    ds_match = save_ds(ds_match, dst_path, encoding = "initiate", label = label)

    return ds_match

def reproject(src_ds, example_ds, dst_path, spatial_interp = "nearest", 
                max_bytes = 2e9, min_times = 10, stack_dim = "time"):

    test_ds = [src_ds, example_ds][np.argmax([src_ds.nbytes, example_ds.nbytes])]

    reproj = choose_reprojecter(test_ds, max_bytes = max_bytes, min_times = min_times, stack_dim = stack_dim)

    if "source" in src_ds.encoding.keys():
        log.info(f"--> Selected `{reproj.__name__}` for reprojection of {os.path.split(src_ds.encoding['source'])[-1]}.").add()
    else:
        log.info(f"--> Selected `{reproj.__name__}` for reprojection.").add()

    ds = reproj(src_ds, example_ds, dst_path, spatial_interp = spatial_interp, stack_dim = stack_dim)
    
    log.sub()

    return ds

def reproject_chunk(src_ds, example_ds, dst_path, spatial_interp = "nearest", stack_dim = "time"):

    das = list()
    variables = dict()
    ncs = list()

    if has_changed(src_ds) or not has_geotransform(src_ds):
        src_path = src_ds.encoding.get("source", dst_path).replace(".nc", "_fixed.nc")
        new_src_ds = src_ds.sortby("y", ascending = False)
        new_src_ds = new_src_ds.rio.write_transform(new_src_ds.rio.transform(recalc=True))
        src_ds = save_ds(new_src_ds, src_path, encoding = "initiate", label = "Correcting src_ds.")
        ncs.append(src_path)
    else:
        src_path = src_ds.encoding["source"]

    resampling = {
                    'nearest': 0,
                    'bilinear': 1,
                    'cubic': 2,
                    'cubic_spline': 3,
                    'lanczos': 4,
                    'average': 5,
                    'mode': 6
                }[spatial_interp]

    example_ds = example_ds.sortby("y", ascending=False)

    vrt_options = {
        'resampling': resampling,
        'crs': example_ds.rio.crs,
        'transform': example_ds.rio.transform(recalc=True),
        'height': example_ds.y.size,
        'width': example_ds.x.size,
        'dtype': "float64",
    }

    for var in src_ds.data_vars:

        part_path = dst_path.replace(".nc", f"_{var}_temp.nc")

        @performance_check
        def _save_warped_vrt(src_path, var, vrt_options, part_path):
            with rasterio.open(f'NETCDF:{src_path}:{var}') as src:
                with WarpedVRT(src, **vrt_options) as vrt:
                    rio_shutil.copy(vrt, part_path, driver='netcdf', creation_options = {"COMPRESS": "DEFLATE"})

        _save_warped_vrt(src_path, var, vrt_options, part_path, label = "Warping VRT to netCDF.")

        ds_part = open_ds(part_path, decode_times = False)

        if stack_dim in src_ds.coords:
            da = ds_part.to_array(stack_dim, name = var).assign_coords({stack_dim: src_ds[stack_dim]})
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

    label = f"Saving reprojected data from {os.path.split(src_ds.encoding['source'])[-1]}:{list(src_ds.data_vars)[0]} ({spatial_interp})."
    ds = save_ds(ds, dst_path, encoding = "initiate", label = label)

    for nc in ncs:
        remove_ds(nc)

    return ds

# if __name__ == "__main__":

#     src_ds = open_ds(r"/Users/hmcoerver/Local/test_data/SENTINEL3/SL_2_LST___.nc")
#     example_ds = open_ds(r"/Users/hmcoerver/Local/test_data/SENTINEL2/S2MSI2A.nc")
#     dst_path = r"/Users/hmcoerver/Local/test_data/output_test.nc"
#     spatial_interp = "nearest"
#     var = "lst"

#     ds = reproject_chunk(src_ds, example_ds, dst_path, spatial_interp = "nearest")
