"""
Functions to prepare input for `pywapor.et_look`, more specifically to
group various parameters in time to create composites. 
"""
import os
import numpy as np
import xarray as xr
import pandas as pd
from pywapor.general.logger import log
import numpy as np
from pywapor.general.processing_functions import save_ds, open_ds
import os
from pywapor.general.reproject import align_pixels
from pywapor.enhancers.apply_enhancers import apply_enhancer
import xarray as xr
import pandas as pd
import dask
import types
import functools

dask.config.set(**{'array.slicing.split_large_chunks': True})

def add_times(ds, bins, composite_type):
    """Add times to the time coordinates, so that every bin has at least one
    datapoint.

    Parameters
    ----------
    ds : xr.Dataset
        Datasat for which to check empty bins.
    bins : list
        List of np.datetime64's which are the boundaries of the groups into
        which the variables will grouped.
    composite_type : {"min" | "max" | "mean"}
        Type of composites that will be created based on the data inside each bin.

    Returns
    -------
    xr.Dataset
        Dataset to which time coordinates have been added to assure no empty bins
        exist.
    """
    
    try:
        bin_count = ds.time.groupby_bins("time", bins).count()
        empty_bins = bin_count.sel({"time_bins": bin_count.isnull()}).time_bins
    except ValueError as e:
        if "None of the data falls within bins" in str(e):
            x = [pd.Interval(pd.Timestamp(x0), 
                            pd.Timestamp(x1)) for x0, x1 in zip(bins[:-1], bins[1:])]
            empty_bins = xr.DataArray(x, {"time_bins": x})
        else:
            raise e

    if composite_type == "mean":
        new_t = [x.mid for x in empty_bins.values]
    else:
        new_t1 = [x.left + pd.Timedelta(seconds=1) for x in empty_bins.values]
        new_t2 = [x.right - pd.Timedelta(seconds=1) for x in empty_bins.values]
        new_t = new_t1 + new_t2

    if len(new_t) > 0:
        ds = xr.merge([ds, xr.Dataset({"time": new_t})]).sortby("time")

    return ds

def time_bins(timelim, bin_length):
    """Based on the time limits and the bin length, create the bin boundaries.

    Parameters
    ----------
    timelim : list
        Period for which to prepare data.
    bin_length : {int | "DEKAD"}
        Length of the bins in days or "DEKAD" for dekadal bins.

    Returns
    -------
    list
        List of np.datetime64's which are the boundaries of the groups into
        which the variables will grouped.
    """
    sdate = timelim[0]
    edate = timelim[1]
    if bin_length == "DEKAD":
        dekad1 = pd.date_range(sdate - pd.to_timedelta(35, unit='d'), edate + pd.to_timedelta(35, unit='d'), freq = "MS")
        dekad2 = dekad1 + pd.Timedelta(10, unit='D')
        dekad3 = dekad1 + pd.Timedelta(20, unit='D')
        dates = np.sort(np.array([dekad1, dekad2, dekad3]).flatten())
        big_interval = pd.Interval(pd.Timestamp(sdate), pd.Timestamp(edate))
        intervals = [pd.Interval(pd.Timestamp(x), pd.Timestamp(y)) for x,y in zip(dates[:-1], dates[1:])]
        out = np.unique([[np.datetime64(x.right, "ns"), np.datetime64(x.left, "ns")] for x in intervals if x.overlaps(big_interval)]).flatten()
    else:
        days = (edate - sdate).days
        no_periods = int(days // bin_length + 1)
        dates = pd.date_range(sdate, periods = no_periods + 1 , freq = f"{bin_length}D")
        out = dates.to_numpy()
    return out

def diags(diagnostics, dss, var):
    dss_diag = list()
    for i, ds in enumerate(dss):
        if isinstance(ds, str):
            ds = open_ds(ds)
        xs = [v[1] for v in diagnostics.values()]
        ys = [v[0] for v in diagnostics.values()]
        ds = ds.sel(x = xs, y = ys, method = "nearest")
        ds[f"{var}_source"] = xr.where(ds[var].notnull(), i, -1)
        dss_diag.append(ds)
    return dss_diag

def create_diags_attrs(srcs):
    attr_dict = dict()
    for i, v in enumerate(srcs):
        if isinstance(v[0], str):
            attr_dict[str(i)] = ".".join(v)
        elif isinstance(v[0], functools.partial):
            attr_dict[str(i)] = v[0].func.__name__
        elif isinstance(v[0], types.FunctionType):
            attr_dict[str(i)] = v[0].__name__
    return attr_dict

def main(dss, sources, example_source, bins, folder, enhancers, cleanup = True):
    """Create composites for variables contained in the 'xr.Dataset's in 'dss'.

    Parameters
    ----------
    dss : dict
        Keys are tuples of ('source', 'product_name'), values are xr.Dataset's 
        which will be aligned along the time dimensions.
    sources : dict
        Configuration for each variable and source.
    example_source : tuple, optional
        Which source to use for spatial alignment, overrides product selected
        through sources, by default None.
    bins : list
        List of 'np.datetime64's which are the boundaries of the groups into
        which the variables will grouped.
    folder : str
        Path to folder in which to store (intermediate) data.
    enhancers : list | "default", optional
        Functions to apply to the xr.Dataset before creating the final
        output, by default "default".
    diagnostics : dict, optional
        Dictionary with coordinates and point-labels for which graphs can be 
        created.

    Returns
    -------
    xr.Dataset
        Dataset with variables grouped into composites.
    """
    # Open unopened netcdf files.
    dss = {**{k: open_ds(v) for k, v in dss.items() if isinstance(v, str)}, 
            **{k:v for k,v in dss.items() if not isinstance(v, str)}}

    final_path = os.path.join(folder, "et_look_in.nc")

    dss2 = list()
    temp_files2 = list()

    compositers = {
        "mean": xr.DataArray.mean,
        "min": xr.DataArray.min,
        "max": xr.DataArray.max,
    }

    for var, config in sources.items():

        spatial_interp = config["spatial_interp"]
        temporal_interp = config["temporal_interp"]
        composite_type = config["composite_type"]
        
        # Align pixels of different products for a single variable together.
        dss_part = [ds[[var]] for ds in dss.values() if var in ds.data_vars]
        dss1, temp_files1 = align_pixels(dss_part, folder, spatial_interp, fn_append = "_step1")

        # Combine different source_products (along time dimension).
        if np.all(["time" in x[var].dims for x in dss1]):
            ds = xr.combine_nested(dss1, concat_dim = "time").chunk({"time": -1}).sortby("time").squeeze()
        elif np.all(["time" not in x[var].dims for x in dss1]):
            ds = xr.concat(dss1, "stacked").median("stacked")
            if len(dss1) > 1:
                log.warning(f"--> Multiple ({len(dss1)}) sources for time-invariant `{var}` found, reducing those with 'median'.")
        else:
            ds = xr.combine_nested([x for x in dss1 if "time" in x[var].dims], concat_dim = "time").chunk({"time": -1}).sortby("time").squeeze()
            log.warning(f"--> Both time-dependent and time-invariant data found for `{var}`, dropping time-invariant data.")

        if "time" in ds.dims:

            if temporal_interp:

                # When combining different products, it is possible to have images 
                # on the exact same date & time. In that case, the median of those images
                # is used. So gaps in one image are filled in with data from the other
                # image(s).
                if np.unique(ds.time).size != ds.time.size:
                    groups = ds.groupby(ds["time"])
                    ds = groups.median(dim = "time")
                    ds = ds.chunk({"time": -1, "y": "auto", "x": "auto"})
                    log.warning(f"--> Multiple `{var}` images for an identical datetime found, reducing those with 'median'.")

                ds = add_times(ds, bins, composite_type = composite_type)
                ds = ds.interpolate_na(dim="time", method = temporal_interp)

            # Make composites.
            ds[var] = compositers[composite_type](ds[var].groupby_bins("time", bins, labels = bins[:-1]))

        # Save output
        dst_path = os.path.join(folder, f"{var}_bin.nc")
        ds = save_ds(ds, dst_path, label = f"Compositing `{var}` ({composite_type}).")
        dss2.append(ds)
        temp_files2.append(ds.encoding["source"])

        for nc in temp_files1:
            if cleanup:
                os.remove(nc)

    # Align all the variables together.
    example_ds = dss[example_source]
    spatial_interps = [sources[list(x.data_vars)[0]]["spatial_interp"] for x in dss2]
    dss3, temp_files3 = align_pixels(dss2, folder, spatial_interps, example_ds, stack_dim = "time_bins", fn_append = "_step2")

    for nc in temp_files2:
        if cleanup:
            os.remove(nc)
    
    ds = xr.merge(dss3)

    # Apply product specific functions.
    for func in enhancers:
        ds, label = apply_enhancer(ds, var, func)
        log.info(label)

    while os.path.isfile(final_path):
        final_path = final_path.replace(".nc", "_.nc")

    ds = save_ds(ds, final_path, encoding = "initiate", label = f"Creating merged file `{os.path.split(final_path)[-1]}`.")

    for nc in temp_files3:
        if cleanup:
            os.remove(nc)

    return ds

if __name__ == "__main__":

    cleanup = False

#     import datetime

#     sources = {
#         "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
#         "r0":           [("MODIS", "MCD43A3.061")],
#         "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
#         "z":            [("SRTM", "30M")],
#         "p":            [("CHIRPS", "P05")],
#         "ra":           [("MERRA2", "M2T1NXRAD.5.12.4")],
#         "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
#         # "t_air_max":    [("MERRA2", "M2I1NXASM.5.12.4")],
#     }

#     diagnostics = { # label          # lat      # lon
#                     "water":	    (29.44977,	30.58215),
#                     "desert":	    (29.12343,	30.51222),
#                     "agriculture":	(29.32301,	30.77599),
#                     "urban":	    (29.30962,	30.84109),
#                     }

#     folder = r"/Users/hmcoerver/Downloads/pywapor_test"
#     latlim = [28.9, 29.7]
#     lonlim = [30.2, 31.2]
#     timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
#     example_source = ("MODIS", "MOD13Q1.061")
#     bin_length = 4

#     bins = time_bins(timelim, bin_length)
#     dss = collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

#     ds = main(dss, sources, example_source, bins, folder, 
#                 diagnostics = None)

#     if diagnostics:
#         _ = main(dss, sources, example_source, bins, folder, diagnostics = diagnostics)