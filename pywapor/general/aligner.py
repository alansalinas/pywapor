"""
Functions to prepare input for `pywapor.se_root`, more specifically to
interpolate various parameters in time to match with land-surface-temperature
times. 
"""

from pywapor.general.processing_functions import save_ds, open_ds
from pywapor.general.reproject import reproject
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.logger import log
import os
import numpy as np
import xarray as xr
from rasterio import CRS

def is_aligned(ds, example_ds):
    return ds.equals(example_ds)

def get_pixel_sizes(dss):
    # Check CRSs of datasets.
    crss = [v.rio.crs.to_epsg() for v in dss]
    # Count occurence of the different CRSs.
    uniqs, counts = np.unique(crss, return_counts=True)
    # Pick the dominant CRS.
    crs = uniqs[np.argmax(counts)]
    # Reproject to common CRS.
    dss = [ds.rio.reproject(CRS.from_epsg(crs)) for ds in dss]
    return [np.abs(np.prod(v.rio.resolution())) for k, v in dss]

def align_pixels(dss, folder, spatial_interp = "nearest", example_ds = None, fn_append = ""):

    temp_files = list()

    if isinstance(spatial_interp, str):
        spatial_interp = [spatial_interp] * len(dss)
    assert len(dss) == len(spatial_interp)

    if len(dss) == 1 and isinstance(example_ds, type(None)):
        dss1 = [dss[0]]
    else:
        if isinstance(example_ds, type(None)):
            example_ds = dss[np.argmin(get_pixel_sizes(dss))]
        dss1 = list()
        for i, (spat_interp, ds_part) in enumerate(zip(spatial_interp, dss)):
            if not ds_part.equals(example_ds):
                var_str = "_".join(ds_part.data_vars)
                dst_path = os.path.join(folder, f"{var_str}_x{i}{fn_append}.nc")
                ds_part = reproject(ds_part, example_ds, dst_path, spatial_interp = spat_interp)
                temp_files.append(ds_part.encoding["source"])
            dss1.append(ds_part)

    return dss1, temp_files

def main(dss, sources, example_source, folder, enhancers, example_t_vars = ["lst"]):
    """Aligns the datetimes in de `dss` xr.Datasets with the datetimes of the 
    dataset with variable `example_t_var`.

    Parameters
    ----------
    dss : dict
        Keys are tuples of (`source`, `product_name`), values are xr.Dataset's 
        or paths (str) to netcdf files, which will be aligned along the time dimensions.
    sources : dict
        Configuration for each variable and source.
    example_source : tuple, optional
        Which source to use for spatial alignment, overrides product selected
        through sources, by default None.
    folder : str
        Path to folder in which to store (intermediate) data.
    enhancers : list | "default", optional
        Functions to apply to the xr.Dataset before creating the final
        output, by default "default".
    example_t_var : str, optional
        Which variable to align the other datasets to in the time dimension, by default "lst".

    Returns
    -------
    xr.Dataset
        Dataset in which all variables have been interpolated to the same times.
    """
    # Open unopened netcdf files.
    dss = {**{k: open_ds(v) for k, v in dss.items() if isinstance(v, str)}, 
            **{k:v for k,v in dss.items() if not isinstance(v, str)}}

    # Determine final output path.
    final_path = os.path.join(folder, "se_root_in.nc")

    # Make list to store intermediate datasets.
    dss2 = list()

    # Make inventory of all variables.
    variables = np.unique([list(ds.data_vars) for ds in dss.values()]).tolist()
    variables = example_t_vars + [x for x in variables if x not in example_t_vars]

    # Create variable to store times to interpolate to.
    example_time = None

    # Loop over the variables
    for var in variables:

        config = sources[var]
        spatial_interp = config["spatial_interp"]
        temporal_interp = config["temporal_interp"]
        # srcs = [(x["source"], x["product_name"]) for x in config["products"]]

        # Align pixels of different products for a single variable together.
        dss_part = [ds[[var]] for ds in dss.values() if var in ds.data_vars]
        dss1, temp_files1 = align_pixels(dss_part, folder, spatial_interp, fn_append = "_step1")

        # Combine different source_products (along time dimension).
        ds = xr.combine_nested(dss1, concat_dim = "time").chunk({"time": -1}).sortby("time").squeeze()

        # TODO fix this in VIIRSL1
        ds = ds.drop_duplicates("time", "last").dropna("time", how = "all")

        temp_files2 = list()
        if var in example_t_vars:
            if isinstance(example_time, type(None)):
                example_time = ds["time"]
            else:
                example_time = xr.concat([example_time, ds["time"]], dim = "time").drop_duplicates("time").sortby("time")
            dss2.append(ds)
        elif "time" in ds[var].dims:
            lbl = f"Aligning times in `{var}` ({ds.time.size}) with `{'` and `'.join(example_t_vars)}` ({example_time.time.size}, {temporal_interp})."
            ds = ds.interpolate_na(dim = "time", method = temporal_interp).ffill("time").bfill("time")
            ds = ds.interp_like(example_time, method = temporal_interp)
            dst_path = os.path.join(folder, f"{var}_i.nc")
            ds = save_ds(ds, dst_path, encoding = "initiate", label = lbl)
            dss2.append(ds)
            temp_files2.append(ds.encoding["source"])

        for nc in temp_files1:
            os.remove(nc)

    # Align all the variables together.
    example_ds = dss[example_source]
    spatial_interps = [sources[list(x.data_vars)[0]]["spatial_interp"] for x in dss2]
    dss3, temp_files3 = align_pixels(dss2, folder, spatial_interps, example_ds, fn_append = "_step2")

    for nc in temp_files2:
        os.remove(nc)

    # Merge everything.
    ds = xr.merge(dss3)

    # Apply product specific functions.
    for func in enhancers:
        ds, label = apply_enhancer(ds, var, func)
        log.info(label)

    if os.path.isfile(final_path):
        final_path = final_path.replace(".nc", "_.nc")

    encoding = {v: {"zlib": True} for v in ds.data_vars}
    ds = save_ds(ds, final_path, encoding = encoding,
                    label = f"Creating merged file `{os.path.split(final_path)[-1]}`.")

    for nc in temp_files3:
        os.remove(nc)

    return ds

if __name__ == "__main__":

    import numpy as np
    import pandas as pd
    from pywapor.general.processing_functions import create_dummy_ds

    # dss = {
    #     ("SENTINEL2", "S2MSI2A"): r"/Users/hmcoerver/Local/test_data/SENTINEL2/S2MSI2A.nc",
    #     ("SENTINEL3", "SL_2_LST___"): r"/Users/hmcoerver/Local/test_data/SENTINEL3/SL_2_LST___.nc",
    #     ("VIIRSL1", "VNP02IMG"): r"/Users/hmcoerver/Local/test_data/VIIRSL1/VNP02IMG.nc",
    # }
    # example_source = ("SENTINEL2", "S2MSI2A")

    sources = {
            "ndvi": {"spatial_interp": "nearest", "temporal_interp": "linear"},
            "lst": {"spatial_interp": "nearest", "temporal_interp": "linear"},
            "bt": {"spatial_interp": "nearest", "temporal_interp": "linear"},
                }

    example_source = ("source1", "productX")
    folder = r"/Users/hmcoerver/Local/aligner_test"
    enhancers = list()
    example_t_vars = ["lst", "bt"]

    chunkers = [(1, 500, 500), (50, 500, 500), (50, 1000, 1000), (50, 3000, 3000)]

    from pywapor.general.logger import log, adjust_logger
    adjust_logger(True, folder, "INFO")

    for i, chunks in enumerate(chunkers):

        dss = {
            ("source1", "productX"): create_dummy_ds(["ndvi"], shape = (10, 2500, 2500), sdate = "2022-02-02", edate = "2022-02-13", fp = os.path.join(folder, f"ndvi_in_test_{i}.nc"), chunks = chunks),
            ("source2", "productY"): create_dummy_ds(["lst"], shape = (23, 500, 500), sdate = "2022-02-01", edate = "2022-02-9", fp = os.path.join(folder, f"lst_in_test_{i}.nc"), min_max = [280, 320], precision=0, chunks = chunks),
            ("source2", "productZ"): create_dummy_ds(["bt"], shape = (14, 420, 420), sdate = "2022-02-04", edate = "2022-02-15", fp = os.path.join(folder, f"bt_in_test_{i}.nc"), min_max = [290, 330], precision=0, chunks = chunks),
        }

        log.info("-".join([str(x) for x in list(dss.values())[0][list(list(dss.values())[0].data_vars)[0]].encoding["chunksizes"]]))

        ds = main(dss, sources, example_source, folder, enhancers, example_t_vars = example_t_vars)

        os.remove(os.path.join(folder, "se_root_in.nc"))
        # os.rename(os.path.join(folder, "log.txt"), os.path.join(folder, f"log_{'-'.join([str(x) for x in chunks])}.txt"))