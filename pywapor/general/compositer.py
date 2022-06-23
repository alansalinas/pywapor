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
from pywapor.general.reproject import reproject
import os
from pywapor.post_et_look import plot_composite
from pywapor.enhancers.apply_enhancers import apply_enhancer
import xarray as xr
import pandas as pd
import dask
import warnings
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

def main(dss, sources, example_source, bins, folder, enhancers,
                diagnostics = None):
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
    warnings.filterwarnings("ignore", message="All-NaN slice encountered")

    final_path = os.path.join(folder, "et_look_in.nc")

    example_ds = open_ds(dss[example_source], "all")

    dss2 = list()

    compositers = {
        "mean": xr.core.groupby.DataArrayGroupBy.mean,
        "min": xr.core.groupby.DataArrayGroupBy.min,
        "max": xr.core.groupby.DataArrayGroupBy.max,
    }

    for var, config in sources.items():

        spatial_interp = config["spatial_interp"]
        temporal_interp = config["temporal_interp"]
        composite_type = config["composite_type"]
        
        srcs = list()
        for x in config["products"]:
            if isinstance(x["source"], str):
                srcs.append((x["source"], x["product_name"]))
            elif isinstance(x["source"], types.FunctionType):
                srcs.append((x["source"].__name__, x["product_name"]))
            elif isinstance(x["source"], functools.partial):
                srcs.append((x["source"].func.__name__, x["product_name"]))

        # Align pixels.
        dst_path = os.path.join(folder, f"{var}.nc")
        if diagnostics:
            dst_path = os.path.join(folder, f"{var}_diags.nc")

        if os.path.isfile(dst_path):
            dss2.append(open_ds(dst_path, "all"))
            continue

        dss1 = [reproject(open_ds(dss[src])[[var]], example_ds, dst_path.replace(".nc", f"_x{i}.nc"), spatial_interp = spatial_interp) for i, src in enumerate(srcs)]

        # TODO FIX this at collect-level (i.e. remove weirds coords using `.drop_vars`)! --> MODIS_Grid_16DAY_250m_500m_VI_eos_cf_projection
        dss1 = [x.drop_vars([y for y in list(x.coords) if y not in ["x", "y", "time", "spatial_ref"]]) for x in dss1]

        if diagnostics:
            dss1 = diags(diagnostics, dss1, var)

        # Combine different source_products (along time dimension).
        ds = xr.combine_nested(dss1, concat_dim = "time").sortby("time")

        if ds.time.size == 1:
            ds = ds.squeeze("time")

        if diagnostics:
            ds[f"{var}_source"].attrs = create_diags_attrs(srcs)
            ds.attrs = {"bin_end": str(bins[-1]), "comp_type": str(composite_type)}

        if "time" in ds.dims:

            ds = ds.chunk({"time": -1, "y": "auto", "x": "auto"})

            if temporal_interp:
                ds = add_times(ds, bins, composite_type = composite_type)
            
            if diagnostics:
                ds[f"{var}_source"] = ds[f"{var}_source"].fillna(255)
                ds[f"{var}_values"] = ds[var]
                log.info(f"--> Compositing `{var}` ({composite_type}) (diagnostic).")
            else:
                log.info(f"--> Compositing `{var}` ({composite_type}).")

            if temporal_interp:

                # When combining different products, it is possible to have images 
                # on the exact same date & time. In that case, the median of those images
                # is used. So gaps in one image are filled in with data from the other
                # image(s).
                if np.unique(ds.time).size != ds.time.size:
                    groups = ds.groupby(ds["time"])
                    ds = groups.median(dim = "time")
                    ds = ds.chunk({"time": -1, "y": "auto", "x": "auto"})
                    log.warning(f"--> Multiple `{var}` images for the same date & time found, reducing those with 'median'.")

                ds = ds.interpolate_na(dim="time", method = temporal_interp)

            # Make composites.
            ds[var] = compositers[composite_type](ds[var].groupby_bins("time", bins, labels = bins[:-1]))

            if diagnostics:
                log.info(f"--> Creating graph for `{var}`.")
                plot_composite(ds, diagnostics, out_folder = os.path.join(folder, "GRAPHS"))

        # Save output
        dss2.append(save_ds(ds, dst_path, decode_coords = "all"))

        for nc in dss1:
            os.remove(nc.encoding["source"])

    if diagnostics:
        for nc in dss2:
            os.remove(nc.encoding["source"])
        return None

    ds = xr.merge(dss2, compat = "override")

    # Apply product specific functions.
    for func in enhancers:
        ds, label = apply_enhancer(ds, var, func)
        log.info(label)

    ds = ds.drop_vars([x for x in ds.coords if (x not in ds.dims) and (x != "spatial_ref")])
    if "time" in list(ds.variables):
        ds = ds.drop_vars("time")
    
    while os.path.isfile(final_path):
        final_path = final_path.replace(".nc", "_.nc")
    ds = save_ds(ds, final_path, decode_coords = "all")

    for nc in dss2:
        os.remove(nc.encoding["source"])

    return ds

# if __name__ == "__main__":

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