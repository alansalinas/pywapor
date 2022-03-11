#%%
import os
import numpy as np
import xarray as xr
from datetime import datetime as dat
import pandas as pd
from dask.diagnostics import ProgressBar
import pywapor.post_et_look as post_et_look
import pywapor.general.processing_functions as PF
import pywapor.general.pre_defaults as defaults
from pywapor.general.logger import log
from pywapor.enhancers.apply_enhancers import apply_enhancer

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
        lean_output = True, diagnostics = {}, extra_source_enhancements = {}):
    """Create temporal composites from spatio-temporal data.

    Parameters
    ----------
    cmeta : dict
        Gives configuration option for the composite creation. Should contain at least 
        `composite_type`, `temporal_interp`, `spatial_interp` and `var_name` as keys.
    dbs : list
        list of lists, where each sublists contains the paths to geotiff files stemming
        from a single source.
    epochs_info : tuple
        Contains three lists with the index-number, starttime and endtime of each
        epoch.
    temp_folder : str, optional
        Path to folder in which temporary files are stores, providing `None` will
        force the function to do all calculations in memory, by default None.
    example_ds : xr.Dataset
        Dataset used to interpolate data (spatially), by default None.
    lean_output : bool, optional
        Whether or not to remove intermediate variables from final output, 
        used for debugging, by default True.
    diagnostics : dict, optional
        Dictionary with (`lat`, `lon`) tuples as values and (short) point descriptions/names
        as values. Several graphs are created for data at these points-of-interest, by default {}.
    extra_source_enhancements : dict, optional
        Dictionary with extra functions to apply to the data before the composites are created, 
        by default {}

    Returns
    -------
    xr.Dataset
        Dataset with the calculated composites.
    """

    source_enhancements = defaults.source_enhancements_defaults()
    source_enhancements = {**source_enhancements, **extra_source_enhancements}

    ## STEP 0: PREPARATIONS.
    if not isinstance(temp_folder, type(None)):
        fh = dbs[0][0]
        path = os.path.normpath(fh)
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)

    if diagnostics:
        diags = {k:v for k,v in diagnostics.items()}
        if "folder" in diags.keys():
            graph_folder = diags.pop("folder")
        else:
            graph_folder = None
        labels = [None, None, None]
        log.info("--> Calculating diagnostics.")
    else:
        diags = None
        labels = ["--> Resampling datasets.",
                  "--> Interpolating {x:,} of {y:,} missing pixels.",
                  "--> Calculating composites."]

    composite_type = cmeta["composite_type"]
    temporal_interp = cmeta["temporal_interp"]
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

    # # Check if data is static
    if np.all([isinstance(composite_type, type(False)), 
                isinstance(temporal_interp, type(False)),
                len(dbs) == 1,
                len(dbs[0]) == 1]):
        ds = xr.open_dataset(dbs[0][0], engine="rasterio")
        ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
        ds = ds.swap_dims({"x": f"lon", "y": f"lat"})
        ds = ds.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
        ds, fh_intermediate = calculate_ds(ds, None, labels[0], 
                                        encoding = {"sources": {"dtype": "uint8"}})
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

    ## STEP 1: COMBINE DATA AND APPLY ENHANCERS.
    ds = None
    dss = list()
    dbs_names = list()

    # Open tif-files and apply spatial interpolation
    for db in dbs:
        check_geots(db)
        sub_ds = xr.open_mfdataset(sorted(db), concat_dim = "time", engine="rasterio", combine = "nested",
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
        ds[var].attrs = {"sources": dbs_names}

    if not isinstance(temp_folder, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_ts.nc")
    else:
        fh = None

    if temporal_interp:
        chunks = {"lat": "auto", "lon": "auto"}
    else:
        chunks = "auto"

    ds, _ = calculate_ds(ds, fh, labels[0], chunks = chunks, 
                            encoding = {"sources": {"dtype": "uint8"}})

    ## STEP 2: TEMPORAL INTERPOLATION.
    if isinstance(epochs_info, type(None)):
        return ds
    else:
        epochs, epoch_starts, epoch_ends = epochs_info

    da = xr.DataArray(epochs, coords={"time":epoch_starts})
    ds["epochs"] = da.reindex_like(ds["time"], method = "ffill", tolerance = epoch_ends[-1] - epoch_starts[-1])

    if temporal_interp:
        t_new = req_t(ds["epochs"], epochs_info, composite_type)
        ds, x, y = temporal_interpolation(ds, t_new, temporal_interp)
        if isinstance(labels[1], str):
            labels[1] = labels[1].format(x = x[all_variables[0]].values, y = y[all_variables[0]].values)
        if x[all_variables[0]].values > 0:
            ds, _ = calculate_ds(ds, fh, label = labels[1], chunks = "auto")
        else: 
            ds = ds.chunk("auto")

    da = xr.DataArray(epochs, coords={"time":epoch_starts})
    ds["epochs"] = da.reindex_like(ds["time"], method = "ffill", tolerance = epoch_ends[-1] - epoch_starts[-1])

    ## STEP 3: COMPOSITE
    for var in all_variables:
        if composite_type == "max":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).max(skipna = True).rename({"epochs": "epoch"})
        elif composite_type == "mean":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).mean(skipna = True).rename({"epochs": "epoch"})
        elif composite_type == "min":
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).min(skipna = True).rename({"epochs": "epoch"})
        elif isinstance(composite_type, float):
            ds[f"{var}_composite"] = ds[var].groupby(ds["epochs"]).quantile(composite_type, skipna = True).rename({"epochs": "epoch"})
        else:
            log.warning("No valid composite_type selected.")
        ds[f"{var}_composite"].attrs = ds[var].attrs

    ds = ds.reindex({"epoch": epochs_info[0]})

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

    ds, _ = calculate_ds(ds, fh, label = labels[2], chunks = "auto")

    if isinstance(diags, dict):
        for var in all_variables:
            post_et_look.plot_composite(ds, diags, graph_folder, band_name = var)

    # if not isinstance(fh_intermediate, type(None)):
    #     os.remove(fh_intermediate)

    return ds

def temporal_interpolation(ds, t_new, temporal_interp):
    """Fill missing and required values using temporal interpolation. In case the
    `ds` contains only a single point in time, any values at `t_new` are extrapolated
    from that single point.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with (missing) data.
    t_new : np.ndarray
        List of extra points in time to add (and interpolate) to `ds`.
    temporal_interp : {"nearest" | "linear"}
        Type of temporal interpolation to apply.

    Returns
    -------
    xr.Dataset
        Dataset with the resampled and interpolated data.
    xr.Dataset
        Dataset with data on how many points can be interpolated (per variable).
    xr.Dataset
        Dataset with data on how many pixels are missing in total (per variable).

    """
    ds_resampled = xr.merge([ds, xr.Dataset({"time": t_new})])
    variables = list(ds_resampled.variables)
    if "sources" in variables:
        ds_resampled["sources"] = ds_resampled.sources.fillna(0.0)
        variables.remove("sources")
    y = ds_resampled.isnull()
    x = (ds_resampled.ffill("time").notnull() * ds_resampled.bfill("time").notnull() * y).sum()
    if ds.time.size > 1:
        for var in variables:
            if "lon" in ds[var].coords and "lat" in ds[var].coords:
                if x[var].values > 0:
                    ds_resampled[var] = ds_resampled[var].interpolate_na(dim = "time", method = temporal_interp)
    elif ds.time.size == 1:
        ds_resampled = ds_resampled.ffill(dim="time").bfill(dim="time")
    return ds_resampled, x, y.sum()

def calculate_ds(ds, fh = None, label = None, encoding = {}, chunks = None):
    """Operations in xarray are often done lazily and this function forces
    calculation and optionally saves the (intermediate) results to disk.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to do operation on.
    fh : str, optional
        Path to `.nc` file in which to save the calculated dataset. `fh` is never overwritten,
        but instead a "_" appended to the filename if the file already exists.
        When `None`, dataset is calculated but not written to disk, by default None.
    label : str, optional
        Label to log, by default None.
    encoding : dict, optional
        Extra encoding to apply when saving to disk, "zlib" is always applied unless
        explicitly set to False with this variable, by default {}.
    chunks : {str | dict}, optional
        Value to chunk the returned data, by default None.

    Returns
    -------
    xr.Dataset
        Returned and calculated dataset.
    str
        Path to saved `.nc` file.

    """
    if isinstance(fh, str):
        while os.path.isfile(fh): # TODO make overwrite instead, was giving weird results so did this for now.
            fh = fh.replace(".nc", "_.nc")
        if not isinstance(label, type(None)):
            log.info(label)
        folder = os.path.split(fh)[0]
        if not os.path.exists(folder):
            os.makedirs(folder)
        with ProgressBar(minimum = 30):
            varis = list(ds.variables)
            default = {"zlib": True}
            full_encoding = {v: default for v in varis if v not in list(ds.coords)}
            for var in encoding.keys():
                full_encoding[var] = {**encoding[var], **full_encoding[var]}
            ds.to_netcdf(fh, encoding = full_encoding)
        new_ds = xr.open_dataset(fh, chunks = chunks)
    else:
        if not isinstance(label, type(None)):
            log.info(label)
        with ProgressBar(minimum=30):
            new_ds = ds.compute()
    return new_ds, fh

def req_t(ds_epochs, epochs_info, composite_type):
    """Determines which epochs don't have any valid data and calculates at which time
    the dataset should be interpolated to be able to create the epochs composite.

    Parameters
    ----------
    ds_epochs : xr.DataArray
        Array with the epoch numbers for which data does currently exist.
    epochs_info : tuple
        Contains three lists with the index-number, starttime and endtime of each
        epoch.
    composite_type : str
        The type of composite determines how many point will need to be interpolated, e.g.
        to have a `mean`-composite, one point in the middle of the epoch will suffice, while
        for a `max` composite two points in time (at the very start and end of the epoch) need
        to be calculated. 

    Returns
    -------
    np.ndarray
        Array of `pd.Timestamps` at which the dataset needs to be interpolated to be able
        to create all the required composites.
    """
    missing_epochs = [x for x in np.unique(epochs_info[0]) if x not in ds_epochs]
    if composite_type == "mean":
        t_new = epochs_info[1][missing_epochs] + (epochs_info[2][missing_epochs] - epochs_info[1][missing_epochs])/2
    else:
        t1 = epochs_info[1][missing_epochs]
        t2 = epochs_info[2][missing_epochs] - pd.Timedelta(seconds=1)
        t_new = np.concatenate((t1, t2))
    return t_new

def check_geots(files):
    """Check if all `files` have the same geotransform and projection.

    Parameters
    ----------
    files : list
        List of paths to geotiff files that are required to have identical 
        geotransforms and projections.
    """
    # TODO: move to PF?
    ref = PF.get_geoinfo(files[0])
    for fh in files[1:]:
        checker = PF.get_geoinfo(fh)
        assert ref == checker, f"ERROR: {files[0]} does not have same geotransform/projection as {fh}."

def get_units(dss):
    """Get an overview of the units of different variables inside xr.Datasets.

    Parameters
    ----------
    dss : list
        List with xr.Datasets. Each variable in each dataset is checked for a 
        `"unit"` attribute.

    Returns
    -------
    dict
        Dictionary with the units per variable.
    """
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
    """Test whether each variable (as key in `units`) has identical units.

    Parameters
    ----------
    units : dict
        Dictionary with the units per variable, can be generated with `compositer.get_units`.
    strictness : {"low" | "medium" | "high}, optional
        low - Units need to be the same, but 'unknown' units are assumed to be correct.
        med - Units need to be the same and 'unknown' units are assumed to be different.
        high - All units must be known and identical.

    Examples
    --------
    >>> units = {"test1": np.array(["C", "C","unknown"]),
    ...          "test2": np.array(["C", "C", "K"]),
    ...          "test3": np.array(["C", "C", "C"]),
    ...          "test4": np.array(["unknown", "unknown", "unknown"])}
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

#%%

if __name__ == "__main__":

    import pywapor
    import pywapor.general.pre_defaults as defaults
    import glob

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"
    epochs_info = pywapor.pre_et_look.create_dates(startdate, enddate, "DEKAD")
    var = "lulc"
    cmeta = defaults.composite_defaults()[var]

    dbs = [[r"/Users/hmcoerver/pywapor_notebooks/RAW/WAPOR/LULC_WAPOR_-_365-daily_2020.01.01.tif"]]

    temp_folder = os.path.join(project_folder, "temporary")
    example_ds = xr.open_dataset(r"/Users/hmcoerver/pywapor_notebooks/example_ds.nc")
    lean_output = True 
    diagnostics = None

#%% STEP1: COMBINE

    # ds = None
    # dss = list()
    # dbs_names = list()

    # # Open tif-files and apply spatial interpolation
    # for db in dbs:
    #     check_geots(db)
    #     sub_ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", combine = "nested",
    #                                     preprocess = preprocess_func)
    #     sub_ds = sub_ds.rename({"band_data": cmeta["var_name"]})

    #     variable = cmeta["var_name"]
    #     source = sub_ds[variable].source

    #     if isinstance(example_ds, type(None)):
    #         example_ds = sub_ds.isel(time = 0).drop_vars(["time"])
    #     else:
    #         sub_ds = sub_ds.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)

    #     dss.append(sub_ds)

    # for i, sub_ds in enumerate(dss):
    #     if i == 0:
    #         ds = sub_ds
    #     else:
    #         ds = ds.combine_first(sub_ds)

    # with ProgressBar():
    #     ds.to_netcdf(r"/Users/hmcoerver/pywapor_notebooks/step1.nc")

#%%
    # STEP 2: TIME INTERPOLATION

    # import xarray as xr
    # from dask.diagnostics import ProgressBar
    # import time

    # fh = r"/Users/hmcoerver/pywapor_notebooks/step1.nc"
    # ds = xr.open_dataset(fh, chunks = {"lat": "auto", "lon": "auto"})
    
    # new_t = req_t()

    # t1 = time.time_ns() 
    # test2 = ds.interpolate_na(dim="time", method="linear")
    # with ProgressBar():
    #     test2 = test2.compute()
    # t2 = time.time_ns()
    # print(i, (t2-t1) * 10**-9)
    
#%%

    # data = np.random.rand(1,3,3)
    # data[0,...] = np.nan
    # # data[2,2,2] = np.nan
    # # data[-1,...] = np.nan

    # ds = xr.Dataset({"ndvi": (["time", "lat","lon"], data)}, coords = {"time": np.arange(data.shape[0])})

    # x = (ds.ffill("time").notnull() * ds.bfill("time").notnull() * ds.isnull()).sum()

    # ds_interpolated = ds.interpolate_na(dim="time", method="nearest", fill_value = "extrapolate")

    # ds.ffill("time").bfill("time")
# %%
    import xarray as xr
    import os
    import numpy as np
    import glob

    # in_file = r"/Users/hmcoerver/pywapor_notebooks/temporary/lst_ts.nc"

    test_files = glob.glob(r"/Users/hmcoerver/pywapor_notebooks/temporary/*.nc")
    more_test_files = glob.glob(r"/Users/hmcoerver/pywapor_notebooks/level_1/*.nc")

    for in_file in test_files + more_test_files:
        print(os.path.split(in_file)[-1])
        ds = xr.open_dataset(in_file)

        default = {
                    # "zlib": True,
                    # "complevel": 1,
                    "scale_factor": 0.001,
                    "dtype": "int32",
                    "_FillValue": -9999
                }
        
        varis = list(ds.variables)
        encoding = {v: default for v in varis if v not in list(ds.coords)}

        out_file = r"/Users/hmcoerver/pywapor_notebooks/test.nc"
        if os.path.isfile(out_file):
            os.remove(out_file)

        ds.to_netcdf(out_file, encoding = encoding)

        print(f"size is --== {os.path.getsize(out_file)/os.path.getsize(in_file)*100:.1f}% ==-- of original")

        ds_compressed = xr.open_dataset(out_file)
        
        # print("lat_dif", ds_compressed.lat.values[:3] - ds.lat.values[:3])
            
        for var in encoding.keys():

            # print(f"{var}_dif", ds_compressed[var].values[0, :3, 0] - ds[var].values[0, :3, 0])

            if "time" in list(ds[var].coords):
                error = np.sqrt(((ds[var].isel(time=0) - ds_compressed[var].isel(time=0))**2).mean())
            elif "epoch" in list(ds[var].coords):
                error = np.sqrt(((ds[var].isel(epoch=0) - ds_compressed[var].isel(epoch=0))**2).mean())
            
            if var == "sources":
                ...
            else:
                print(error.values)

        ds = ds.close()
        ds_compressed = ds_compressed.close()

        print("")

    
# %%
