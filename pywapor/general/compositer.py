import os
import numpy as np
import xarray as xr
import pandas as pd
from pywapor.general.logger import log
import numpy as np
import pywapor
from pywapor.general.processing_functions import save_ds, open_ds
from pywapor.general.reproject import reproject
from pywapor.collect.downloader import collect_sources
import os
import xarray as xr
import pandas as pd

def add_times(ds, bins, composite_type):
    
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
        new_t1 = [x.left for x in empty_bins.values]
        new_t2 = [x.right - pd.Timedelta(seconds=1) for x in empty_bins.values]
        new_t = new_t1 + new_t2

    if len(new_t) > 0:
        ds = xr.merge([ds, xr.Dataset({"time": new_t})]).sortby("time")

    return ds

def time_bins(timelim, bin_length):
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
        ds[f"{var}_source"] = xr.where(ds[var].notnull(), i, np.nan)
        dss_diag.append(ds)
    return dss_diag

def main(dss, sources, example_source, bins, folder, 
                diagnostics = None):

    example_ds = open_ds(dss[example_source], "all")
    cdefaults = pywapor.general.pre_defaults.composite_defaults()

    dss2 = list()

    compositers = {
        "mean": xr.core.groupby.DataArrayGroupBy.mean,
        "min": xr.core.groupby.DataArrayGroupBy.min,
        "max": xr.core.groupby.DataArrayGroupBy.mean,
    }

    for var, srcs in sources.items():

        spatial_interp = cdefaults[var]["spatial_interp"]
        temporal_interp = cdefaults[var]["temporal_interp"]
        composite_type = cdefaults[var]["composite_type"]

        dst_path = os.path.join(folder, f"{var}.nc")

        # Align pixels.
        dss1 = [reproject(open_ds(dss[src])[[var]], example_ds, dst_path.replace(".nc", f"_{i}.nc"), spatial_interp = spatial_interp) for i, src in enumerate(srcs)]

        if diagnostics:
            dst_path = dst_path.replace(".nc", "_diag.nc")
            dss1 = diags(diagnostics, dss1, var)

        # Combine different source_products (along time dimension).
        ds = xr.combine_nested(dss1, concat_dim = "time").sortby("time").squeeze()

        if ("time" in ds[var].dims):
            ds = ds.chunk({"time": -1})

            # Interpolate (temporal).
            ds = add_times(ds, bins, composite_type = composite_type)
            
            if diagnostics:
                ds[f"{var}_source"] = ds[f"{var}_source"].fillna(255)
                ds[f"{var}_values"] = ds[var]
            
            ds = ds.interpolate_na(dim="time", method = temporal_interp)

            # Make composites.
            ds[var] = compositers[composite_type](ds[var].groupby_bins("time", bins, labels = bins[:-1]))

        # Save output
        dss2.append(save_ds(ds, dst_path, decode_coords = "all"))

        for nc in dss1:
            os.remove(nc.encoding["source"])

    final_path = {
                    False: os.path.join(folder, "et_look_in.nc"),
                    True: os.path.join(folder, "et_look_in_diags.nc")
                }[bool(diagnostics)]

    ds = xr.merge(dss2)
    ds = ds.drop_vars([x for x in ds.coords if (x not in ds.dims) and (x != "spatial_ref")])
    ds = save_ds(ds, final_path, decode_coords = "all")

    for nc in dss2:
        os.remove(nc.encoding["source"])

    return ds

if __name__ == "__main__":

    import datetime

    sources = {
        "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
        "r0":           [("MODIS", "MCD43A3.061")],
        "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
        "z":            [("SRTM", "30M")],
        "p":            [("CHIRPS", "P05")],
        "ra":           [("MERRA2", "M2T1NXRAD.5.12.4")],
        "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        # "t_air_max":    [("MERRA2", "M2I1NXASM.5.12.4")],
    }

    diagnostics = { # label          # lat      # lon
                    "water":	    (29.44977,	30.58215),
                    "desert":	    (29.12343,	30.51222),
                    "agriculture":	(29.32301,	30.77599),
                    "urban":	    (29.30962,	30.84109),
                    }

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
    example_source = ("MODIS", "MOD13Q1.061")
    bin_length = 4

    bins = time_bins(timelim, bin_length)
    dss = collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

    ds = main(dss, sources, example_source, bins, folder, 
                diagnostics = None)

    if diagnostics:
        _ = main(dss, sources, example_source, bins, folder, diagnostics = diagnostics)