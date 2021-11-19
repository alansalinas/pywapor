import os
import numpy as np
import xarray as xr
from datetime import datetime as dat
import pandas as pd
from dask.diagnostics import ProgressBar
import pywapor.post_et_look as pl
import pywapor.general.processing_functions as PF

def main(cmeta, dbs, epochs_info, temp_folder, example_ds = None,
                lean_output = True, diagnostics = None):
    """
    composite_type = "max", "mean", "min"
    temporal_interp = False, "linear", "nearest", "zero", "slinear", "quadratic", "cubic"
    spatial_interp = "nearest", "linear"
    """

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    composite_type = cmeta["composite_type"]
    temporal_interp = cmeta["temporal_interp"]
    spatial_interp = cmeta["spatial_interp"]

    styling = { 1: ("*", "r", 1.0, "MOD13Q1"),
                2: ("o", "g", 1.0, "MYD13Q1"),
                3: ("v", "b", 1.0, "PROBAV"),
                4: ("s", "y", 1.0, "MCD43A3"),
                5: ("*", "purple", 1.0, "CHIRPS.v2.0"),
                6: ("p", "darkblue", 1.0, "GEOS"),
                7: ("h", "gray", 1.0, "MERRA"),
                999: ("P", "orange", 1.0, "-"),
                0: (".", "k", 0.7, "Interp.")}

    # Check if data is static
    if np.all([isinstance(composite_type, type(None)), 
                isinstance(temporal_interp, type(None)),
                len(dbs) == 1,
                len(dbs[0]) == 1]):
        ds = xr.open_dataset(dbs[0][0], engine="rasterio")
        ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
        ds = ds.swap_dims({"x": f"lon", "y": f"lat"})
        ds = ds.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
        ds = ds.rename({"band": "epoch"}).assign_coords({"epoch": [-9999]})
        ds = ds.rename({"band_data": "composite"})
        return ds

    def preprocess_func(ds):

        date_string = ds.encoding["source"].split("_")[-1]
        freq_string = ds.encoding["source"].split("_")[-2]
        if len(date_string.split(".")) > 4:
            date = dat.strptime(date_string, '%Y.%m.%d.%H.%M.tif')
        else:
            date = dat.strptime(date_string, '%Y.%m.%d.tif')
            t_offsets = {"daily": 12, "16-daily": 192, 
                        "daily-min":12, "daily-max":12, "5-daily": 60}
            date = date + pd.Timedelta(hours = t_offsets[freq_string])

        source = ds.encoding["source"].split("_")[-4]

        ds = ds.drop_vars("spatial_ref")
        ds = ds.squeeze("band")
        ds = ds.drop_vars("band")
        ds = ds.expand_dims("time", axis = 0)
        ds = ds.assign_coords({"time": [date]})
        ds = ds.rename_vars({"x": f"lon", "y": f"lat"})
        ds = ds.swap_dims({"x": f"lon", "y": f"lat"})

        sources = {v[3]: k for k, v in styling.items()}

        ds["band_data"] = ds.band_data.expand_dims("source")
        ds = ds.assign_coords({"source": [source]})

        ds["sources"] = np.isfinite(ds.band_data) * sources[source]

        return ds

    ds = None
    dss = list()
    dbs_names = list()

    # Open tif-files and apply spatial interpolation
    for db in dbs:

        check_geots(db)

        sub_ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", combine = "nested",
                                        preprocess = preprocess_func)
        dbs_names.append(sub_ds.source.values[0])

        if isinstance(example_ds, type(None)):
            example_ds = sub_ds.isel(time = 0, 
                                    source = 0).drop_vars(["source", "time"])
        else:
            sub_ds = sub_ds.interp_like(example_ds, method = spatial_interp, kwargs={"fill_value": "extrapolate"},)

        dss.append(sub_ds)

    ds = xr.merge(dss)

    for i, db_name in enumerate(dbs_names):
        if i == 0:
            ds_temp = ds.sel(source = db_name)
        else:
            # Combine datasets, this is where one ds is prioritised over another!
            ds_temp = ds_temp.fillna(ds.sel(source = db_name))

    ds = ds_temp

    fh = os.path.join(temp_folder, "intermediate.nc")

    ds = ds.reindex({"time": ds.time.sortby("time")})
    ds = calculate_ds(ds, fh, "--> Reprojecting datasets.").chunk({"time":-1})

    if temporal_interp:
        ds_resampled = ds.interpolate_na(dim="time")
        ds_resampled = ds_resampled.resample({"time": "12H"}).interpolate(temporal_interp).drop("sources")
        ds_resampled["sources"] = ds.sources
        ds_resampled["sources"] = ds_resampled["sources"].fillna(0.0)
        ds = ds_resampled

    if isinstance(epochs_info, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_ts.nc")
        ds = calculate_ds(ds, fh, "--> Exporting timeseries.")
        return ds

    epochs, epoch_starts, epoch_ends = epochs_info

    da = xr.DataArray(epochs, coords={"time":epoch_starts})
    ds["epochs"] = da.reindex_like(ds["time"], method = "ffill", tolerance = epoch_ends[-1] - epoch_starts[-1])

    if composite_type == "max":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).max().rename({"epochs": "epoch"})
    elif composite_type == "mean":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).mean().rename({"epochs": "epoch"})
    elif composite_type == "min":
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).min().rename({"epochs": "epoch"})
    elif isinstance(composite_type, float):
        ds["composite"] = ds.band_data.groupby(ds["epochs"]).quantile(composite_type).rename({"epochs": "epoch"})
    else:
        print("No valid composite_type selected.")

    ds["epoch_starts"] = xr.DataArray(epoch_starts, coords = {"epoch": epochs})
    ds["epoch_ends"] = xr.DataArray(epoch_ends, coords = {"epoch": epochs})

    checklist = {True: dbs_names + ["Interp."], False: dbs_names}[True] # TODO interp data should not be included when not used
    sources_styling = {str(k): v for k, v in styling.items() if v[3] in checklist}
    ds.attrs = {str(k): str(v) for k, v in {**cmeta, **sources_styling}.items()}

    if diagnostics:
        graph_folder = diagnostics.pop("folder")
        ds_diags = ds.sel(lat = [v[0] for v in diagnostics.values()],
                          lon = [v[1] for v in diagnostics.values()], method = "nearest")
    else:
        ds_diags = None

    if lean_output: # TODO this option can be removed?
        ds = ds.drop_vars(["band_data", "sources", "epochs", "time"])
    else:
        ds = ds.drop_vars(["sources", "epochs"])
   
    fh = os.path.join(temp_folder, f"{cmeta['var_name']}_composite.nc")
    ds = calculate_ds(ds, fh, label = "--> Calculating composites.")

    if not isinstance(diagnostics, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_diags.nc")
        ds_diags = calculate_ds(ds_diags, fh, label = "--> Calculating diagnostics.")
        pl.plot_composite(ds_diags, diagnostics, graph_folder)
        diagnostics["folder"] = graph_folder

    return ds

def calculate_ds(ds, fh, label = None):
    if os.path.isfile(fh):
        os.remove(fh)
    with ProgressBar(minimum = 30):
        if not isinstance(label, type(None)):
            print(label)
        ds.to_netcdf(fh)
    ds.close()
    ds = None
    ds = xr.open_dataset(fh)
    return ds

def check_geots(files):
    ref = PF.get_geoinfo(files[0])
    for fh in files[1:]:
        checker = PF.get_geoinfo(fh)
        assert ref == checker, f"ERROR: {files[0]} does not have same geotransform/projection as {fh}."
