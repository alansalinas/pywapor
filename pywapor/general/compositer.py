import os
import numpy as np
import xarray as xr
from datetime import datetime as dat
import pandas as pd
from dask.diagnostics import ProgressBar
from pywapor.collect.downloader import collect_sources
import pywapor.post_et_look as post_et_look
import pywapor.general.processing_functions as PF
from pywapor.general.logger import log
from pywapor.enhancers.temperature import kelvin_to_celsius
import pywapor.enhancers.lulc as lulc
from functools import partial
from pywapor.enhancers.apply_enhancers import apply_enhancer

def remove_var(ds, var):
    return ds.drop_vars([var])

def preprocess_func(ds):

    date_string = ds.encoding["source"].split("_")[-1]
    freq_string = ds.encoding["source"].split("_")[-2]
    unit_string = ds.encoding["source"].split("_")[-3]
    source_string = ds.encoding["source"].split("_")[-4]

    if len(date_string.split(".")) > 4:
        date = dat.strptime(date_string, '%Y.%m.%d.%H.%M.tif')
    else:
        date = dat.strptime(date_string, '%Y.%m.%d.tif')
        t_offsets = {"daily": 12, "daily-min":12, "daily-max":12} # TODO get rid of daily-min/max
        if freq_string in t_offsets.keys():
            delta = t_offsets[freq_string]
        else:
            delta = float(freq_string.split("-")[0]) * 12
        date = date + pd.Timedelta(hours = delta)

    ds = ds.drop_vars("spatial_ref")
    ds = ds.squeeze("band")
    ds = ds.drop_vars("band")
    ds = ds.expand_dims("time", axis = 0)
    ds = ds.assign_coords({"time": [date]})
    ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
    ds = ds.swap_dims({"x": f"lon", "y": f"lat"})

    ds.band_data.attrs = {"unit": unit_string,
                            "source": source_string}

    return ds

def main(cmeta, dbs, epochs_info, temp_folder = None, example_ds = None,
        lean_output = True, diagnostics = None):
    """
    composite_type = "max", "mean", "min"
    temporal_interp = False, "linear", "nearest", "zero", "slinear", "quadratic", "cubic"
    spatial_interp = "nearest", "linear"
    """

    source_enhancements = {
        ("MERRA2",  "t_air_24"):        [kelvin_to_celsius],
        ("MERRA2",  "t_air_min_24"):    [kelvin_to_celsius],
        ("MERRA2",  "t_air_max_24"):    [kelvin_to_celsius],
        ("GEOS5",   "t_air_24"):        [kelvin_to_celsius],
        ("GEOS5",   "t_air_min_24"):    [kelvin_to_celsius],
        ("GEOS5",   "t_air_max_24"):    [kelvin_to_celsius],
        ("GLOBCOVER", "lulc"):          [partial(lulc.lulc_to_x, out_var = "land_mask", 
                                                    convertor = lulc.globcover_to_land_mask()),
                                        partial(lulc.lulc_to_x, out_var = "rs_min", 
                                                    convertor = lulc.globcover_to_rs_min()),
                                        partial(lulc.lulc_to_x, out_var = "lue_max", 
                                                    convertor = lulc.globcover_to_lue_max()),
                                        partial(lulc.lulc_to_x, out_var = "z_obst_max", 
                                                    convertor = lulc.globcover_to_z_obst_max()),
                                        remove_var,
                                        ],
        ("WAPOR", "lulc"):              [partial(lulc.lulc_to_x, out_var = "land_mask", 
                                                    convertor = lulc.wapor_to_land_mask()),
                                        partial(lulc.lulc_to_x, out_var = "rs_min", 
                                                    convertor = lulc.wapor_to_rs_min()),
                                        partial(lulc.lulc_to_x, out_var = "lue_max", 
                                                    convertor = lulc.wapor_to_lue_max()),
                                        partial(lulc.lulc_to_x, out_var = "z_obst_max", 
                                                    convertor = lulc.wapor_to_z_obst_max()),
                                        remove_var,
                                        ],
    }

    if not isinstance(temp_folder, type(None)):
        fh = dbs[0][0]
        path = os.path.normpath(fh)
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

    if isinstance(diagnostics, dict):
        diags = {k:v for k,v in diagnostics.items()}
        if "folder" in diags.keys():
            graph_folder = diags.pop("folder")
        else:
            graph_folder = None
        labels = [None, None]
        log.info("--> Calculating diagnostics.")
    else:
        diags = None
        labels = ["--> Resampling datasets.",
                    "--> Calculating composites."]

    composite_type = cmeta["composite_type"]
    temporal_interp = cmeta["temporal_interp"]
    if isinstance(cmeta["temporal_interp_freq"], int):
        periods = cmeta["temporal_interp_freq"]
        freq = None
    else:
        periods = None
        freq = cmeta["temporal_interp_freq"]
    spatial_interp = cmeta["spatial_interp"]

    styling = dict()
    markers = ["*", "o", "v", "s", "*", "p", "h"]
    colors =  ["r", "g", "b", "y", "purple", "darkblue", "gray", "orange"]
    for i, fps in enumerate(dbs):
        source = os.path.split(fps[0])[-1].split("_")[1]
        styling[i+1] = (markers[i], colors[i], 1.0, source)
    styling[999] = ("P", "orange", 1.0, "-")
    styling[0] = (".", "k", 0.7, "Interp.")
    sources = {v[3]: k for k, v in styling.items()}

    # Check if data is static
    if np.all([isinstance(composite_type, type(False)), 
                isinstance(temporal_interp, type(False)),
                len(dbs) == 1,
                len(dbs[0]) == 1]):
        ds = xr.open_dataset(dbs[0][0], engine="rasterio")
        ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
        ds = ds.swap_dims({"x": f"lon", "y": f"lat"})
        ds = ds.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
        ds, fh_intermediate = calculate_ds(ds, None, labels[0], 
                                        cast = {"sources": np.uint8})
        # ds = ds.rename({"band": "epoch"}).assign_coords({"epoch": [-9999]})
        ds = ds.isel(band = 0)
        ds = ds.drop_vars(["band"])
        ds = ds.rename({"band_data": f"{cmeta['var_name']}"})

        unit_string = dbs[0][0].split("_")[-3]
        source_string = dbs[0][0].split("_")[-4]

        attributes = {
                        # "unit": unit_string,
                        "sources": [source_string],
                    }

        ds[cmeta['var_name']].attrs = attributes
        return ds

    ds = None
    dss = list()
    dbs_names = list()

    # Open tif-files and apply spatial interpolation
    for db in dbs:
        check_geots(db)
        sub_ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", combine = "nested",
                                        preprocess = preprocess_func)
        sub_ds = sub_ds.rename({"band_data": cmeta["var_name"]})

        variable = cmeta["var_name"]
        source = sub_ds[variable].source

        dbs_names.append(source)
        sub_ds["sources"] = (np.isfinite(sub_ds[variable]) * sources[source])

        if (source, variable) in source_enhancements.keys():
            for enhancer in source_enhancements[(source, variable)]:
                sub_ds, label = apply_enhancer(sub_ds, variable, enhancer, source = source, 
                                log_it = isinstance(diags, type(None)))
                sub_ds, _ = calculate_ds(sub_ds, label = label)

        if isinstance(example_ds, type(None)):
            example_ds = sub_ds.isel(time = 0).drop_vars(["time"])
        else:
            sub_ds = sub_ds.interp_like(example_ds, method = spatial_interp, kwargs={"fill_value": "extrapolate"},)

        if isinstance(diags, dict):
            # Select the pixels that contain the diagnostics points.
            sub_ds = sub_ds.sel(lat = [v[0] for v in diags.values()],
                                lon = [v[1] for v in diags.values()], method = "nearest")
            # Remove pixels that where selected twice, in case diagnostics points
            # lie within the same pixel (usually only necessary for very coarse data).
            sub_ds = sub_ds.isel(lat = np.unique(sub_ds.lat, return_index = True)[1], 
                                 lon = np.unique(sub_ds.lon, return_index = True)[1])

        dss.append(sub_ds)

    results = check_units(get_units(dss), strictness = "med") # TODO: Move to "high"
    assert np.all(results.values()), f"Combining data with different units: {get_units(dss)}"

    # This is where one dataset gets priority over another if there are sources
    # at the exact same time. That's also why we can't simply use xr.merge() here.
    for i, sub_ds in enumerate(dss):
        if i == 0:
            ds = sub_ds
        else:
            ds = ds.combine_first(sub_ds)

    all_variables = list(ds.keys())
    all_variables.remove('sources')

    for var in all_variables:
        ds[var].attrs = {"sources": dbs_names,
                        #  "unit": get_units(dss)[var][0]
                         }

    if not isinstance(temp_folder, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_ts.nc")
    else:
        fh = None

    ds, fh_intermediate = calculate_ds(ds, fh, labels[0], 
                                        cast = {"sources": np.uint8})
    ds = ds.chunk("auto")

    if isinstance(epochs_info, type(None)):
        return ds
    else:
        epochs, epoch_starts, epoch_ends = epochs_info

    if temporal_interp:
        t_new = calc_interpolation_times(epochs_info, freq = freq, periods = periods)
        ds_resampled = xr.merge([ds, xr.Dataset({"time": t_new})])
        ds_resampled["sources"] = ds_resampled.sources.fillna(0.0)
        if ds.time.size > 1:
            ds_resampled = ds_resampled.interpolate_na(dim = "time", method = temporal_interp, 
                                        max_gap = None) # TODO max_gap can be a datetime.timedelta object.
        else:
            ds_resampled = ds_resampled.ffill(dim="time")
            ds_resampled = ds_resampled.bfill(dim="time")
        ds = ds_resampled

    da = xr.DataArray(epochs, coords={"time":epoch_starts})
    ds["epochs"] = da.reindex_like(ds["time"], method = "ffill", tolerance = epoch_ends[-1] - epoch_starts[-1])

    for var in all_variables:
        if composite_type == "max":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).max().rename({"epochs": "epoch"})
        elif composite_type == "mean":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).mean().rename({"epochs": "epoch"})
        elif composite_type == "min":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).min().rename({"epochs": "epoch"})
        elif isinstance(composite_type, float):
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).quantile(composite_type).rename({"epochs": "epoch"})
        else:
            log.warning("No valid composite_type selected.")
        ds[f"{var}_composite"].attrs = ds[var].attrs

    # ds["epoch_starts"] = xr.DataArray(epoch_starts, coords = {"epoch": epochs})
    # ds["epoch_ends"] = xr.DataArray(epoch_ends, coords = {"epoch": epochs})

    ds = ds.assign_coords(epoch_starts = ("epoch", epoch_starts))
    ds = ds.assign_coords(epoch_ends = ("epoch", epoch_ends))

    checklist = {True: dbs_names + ["Interp."], False: dbs_names}[True] # TODO interp data should not be included when not used
    sources_styling = {str(k): v for k, v in styling.items() if v[3] in checklist}
    ds.attrs = {str(k): str(v) for k, v in {**cmeta, **sources_styling}.items()}

    coords_to_keep = ["lon", "lat", "epoch", "epoch_starts", "epoch_ends", "epochs", "time"]
    coords_to_drop = [x for x in list(ds.coords) if x not in coords_to_keep]
    ds = ds.drop_vars(coords_to_drop)

    if lean_output:
        # Keep only the composites.
        ds = ds.drop_vars(all_variables + ["sources", "epochs", "time"])

    if not isinstance(temp_folder, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_composite.nc")
    else:
        fh = None

    ds, fh_composites = calculate_ds(ds, fh, label = labels[1])

    if isinstance(diags, dict):
        for var in all_variables:
            post_et_look.plot_composite(ds, diags, graph_folder, band_name = var)

    # if not isinstance(fh_intermediate, type(None)):
    #     os.remove(fh_intermediate)

    return ds

def calculate_ds(ds, fh = None, label = None, cast = None):
    if isinstance(fh, str):
        while os.path.isfile(fh): # TODO make overwrite instead, was giving weird results so did this for now.
            fh = fh.replace(".nc", "_.nc")
        with ProgressBar(minimum = 30):
            if not isinstance(label, type(None)):
                log.info(label)
            if isinstance(cast, dict):
                for k, v in cast.items():
                    ds[k] = ds[k].astype(v)
            ds.to_netcdf(fh)
        new_ds = xr.open_dataset(fh)
    else:
        with ProgressBar(minimum=30):
            if not isinstance(label, type(None)):
                log.info(label)
            new_ds = ds.compute()
    return new_ds, fh

def check_geots(files):
    ref = PF.get_geoinfo(files[0])
    for fh in files[1:]:
        checker = PF.get_geoinfo(fh)
        assert ref == checker, f"ERROR: {files[0]} does not have same geotransform/projection as {fh}."

def calc_interpolation_times(epochs_info, freq = "2D", periods = None):
    new_t = np.array([], dtype = np.datetime64)
    for epoch_start, epoch_end in zip(epochs_info[1], epochs_info[2]):
        if isinstance(periods, type(None)):
            part_freq = freq
        else:
            delta = epoch_end - epoch_start
            if isinstance(delta, np.timedelta64):
                part_freq = f"{int(delta/periods)}N"
            else:
                part_freq = f"{int(delta.value/periods)}N"
        new_part_t = pd.date_range(epoch_start, epoch_end, freq = part_freq)
        new_t = np.append(new_t, new_part_t)
    new_t = np.sort(np.unique(new_t))
    return new_t

def get_units(dss):
    units = dict()
    for sub_ds in dss:
        variables = list(sub_ds.keys())
        for var in variables:
            if hasattr(sub_ds[var], "unit"):
                unit = sub_ds[var].unit
            else:
                unit = "unknown"
            if var in units.keys():
                units[var] = np.append(units[var], unit)
            else:
                units[var] = np.array([unit])
    return units

def check_units(units, strictness = "low"):
    """[summary]

    Parameters
    ----------
    dss : [type]
        [description]
    strictness : str, optional
        low - Units need to be the same, but 'unknown' units are assumed to be correct.
        med - Units need to be the same and 'unknown' units are assumed to be different.
        high - All units must be known and identical.

    Examples
    --------
    >>> units = {"test1": np.array(["C", "C","unknown"]),
             "test2": np.array(["C", "C", "K"]),
             "test3": np.array(["C", "C", "C"]),
             "test4": np.array(["unknown", "unknown", "unknown"])}

    >>> check_units(units)
    {'test1': True, 'test2': False, 'test3': True, 'test4': True}
    >>> check_units(units, strictness = "med")
    {'test1': False, 'test2': False, 'test3': True, 'test4': True}
    >>> check_units(units, strictness = "high")
    {'test1': False, 'test2': False, 'test3': True, 'test4': False}
    """
    results = dict()
    if strictness == "low":
        for k, v in units.items():
            check = np.unique(v[v!="unknown"]).size <= 1
            results[k] = check
    if strictness == "med":
        for k, v in units.items():
            check = np.unique(v).size == 1
            results[k] = check
    if strictness == "high":
        for k, v in units.items():
            check = np.unique(v).size == 1 and "unknown" not in v
            results[k] = check
    return results

if __name__ == "__main__":

    import pywapor
    import pywapor.general.pre_defaults as defaults

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"

    epochs_info = pywapor.pre_et_look.create_dates(startdate, enddate, "DEKAD")

    var = "ndvi"

    cmeta = defaults.composite_defaults()[var]

    sources = ["MOD13", "MYD13"]

    dl_args = {
                "Dir": os.path.join(project_folder, "RAW"), 
                "latlim": latlim, 
                "lonlim": lonlim, 
                "Startdate": startdate, 
                "Enddate": enddate,
                }

    dbs = collect_sources(var, sources, dl_args)

    temp_folder = os.path.join(project_folder, "temporary")
    example_ds = None
    lean_output = True 
    diagnostics = None

    # ds = main(cmeta, dbs, epochs_info, temp_folder = temp_folder, example_ds = example_ds,
    #     lean_output = lean_output, diagnostics = diagnostics)

    # diagnostics = { # label          # lat      # lon
    #                 "water":	    (29.44977,	30.58215),
    #                 "desert":	    (29.12343,	30.51222),
    #                 "agriculture":	(29.32301,	30.77599),
    #                 "urban":	    (29.30962,	30.84109),
    #                 }

    # temp_folder = None
    # lean_output = False

    # ds_diags = main(cmeta, dbs, epochs_info, temp_folder = temp_folder, example_ds = example_ds,
    #     lean_output = lean_output, diagnostics = diagnostics)

