import os
import numpy as np
import xarray as xr
from datetime import datetime as dat
import pandas as pd
import datetime
from dask.diagnostics import ProgressBar
import pywapor.post_et_look as post_et_look
import pywapor.general.processing_functions as PF
from pywapor.general.logger import log

def main(cmeta, dbs, epochs_info, temp_folder = None, example_ds = None,
        lean_output = True, diagnostics = None):
    """
    composite_type = "max", "mean", "min"
    temporal_interp = False, "linear", "nearest", "zero", "slinear", "quadratic", "cubic"
    spatial_interp = "nearest", "linear"
    """

    if not isinstance(temp_folder, type(None)):
        fh = dbs[0][0]
        path = os.path.normpath(fh)
        parts = path.split(os.sep)
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

        return ds

    ds = None
    dss = list()
    dbs_names = list()

    # Open tif-files and apply spatial interpolation
    for db in dbs:
        check_geots(db)
        sub_ds = xr.open_mfdataset(db, concat_dim = "time", engine="rasterio", combine = "nested",
                                        preprocess = preprocess_func)
        if isinstance(example_ds, type(None)):
            example_ds = sub_ds.isel(time = 0).drop_vars(["time"])
        else:
            sub_ds = sub_ds.interp_like(example_ds, method = spatial_interp, kwargs={"fill_value": "extrapolate"},)

        source = db[0].split("_")[-4]
        dbs_names.append(source)
        sub_ds["sources"] = (np.isfinite(sub_ds.band_data) * sources[source])

        if isinstance(diags, dict):
            sub_ds = sub_ds.sel(lat = [v[0] for v in diags.values()],
                        lon = [v[1] for v in diags.values()], method = "nearest")

        dss.append(sub_ds)

    # This is where one dataset gets priority over another if there are sources
    # at the exact same time. That's also why we can't use xr.merge() here.
    for i, sub_ds in enumerate(dss):
        if i == 0:
            ds = sub_ds
        else:
            ds = ds.combine_first(sub_ds)

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
        t_new = np.sort(np.unique(np.append(ds.time, [epoch_starts, epoch_ends, epoch_ends - datetime.timedelta(seconds=1)])))
        ds_resampled = xr.merge([ds, xr.Dataset({"time": t_new})])
        ds_resampled["sources"] = ds_resampled.sources.fillna(0.0)
        ds_resampled = ds_resampled.interpolate_na(dim = "time", method = temporal_interp, 
                                    max_gap = None) # TODO max_gap can be a datetime.timedelta object.
        ds = ds_resampled

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
        log.warning("No valid composite_type selected.")

    ds["epoch_starts"] = xr.DataArray(epoch_starts, coords = {"epoch": epochs})
    ds["epoch_ends"] = xr.DataArray(epoch_ends, coords = {"epoch": epochs})

    checklist = {True: dbs_names + ["Interp."], False: dbs_names}[True] # TODO interp data should not be included when not used
    sources_styling = {str(k): v for k, v in styling.items() if v[3] in checklist}
    ds.attrs = {str(k): str(v) for k, v in {**cmeta, **sources_styling}.items()}

    if lean_output:
        ds = ds.drop_vars(["band_data", "sources", "epochs", "time"])

    if not isinstance(temp_folder, type(None)):
        fh = os.path.join(temp_folder, f"{cmeta['var_name']}_composite.nc")
    else:
        fh = None

    ds, fh_composites = calculate_ds(ds, fh, label = labels[1])

    if isinstance(diags, dict):
        post_et_look.plot_composite(ds, diags, graph_folder)

    if not isinstance(fh_intermediate, type(None)):
        os.remove(fh_intermediate)

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
        ds.close()
        ds = None
        ds = xr.open_dataset(fh)
    else:
        with ProgressBar(minimum=30):
            if not isinstance(label, type(None)):
                log.info(label)
            ds = ds.compute()
    return ds, fh

def check_geots(files):
    ref = PF.get_geoinfo(files[0])
    for fh in files[1:]:
        checker = PF.get_geoinfo(fh)
        assert ref == checker, f"ERROR: {files[0]} does not have same geotransform/projection as {fh}."

# if __name__ == "__main__":

#     import pywapor

#     project_folder = r"/Users/hmcoerver/On My Mac/pyWAPOR/example_data"
#     latlim = [28.9, 29.7]
#     lonlim = [30.2, 31.2]
#     startdate = "2021-06-01"
#     enddate = "2021-08-01"

#     epochs_info = pywapor.pre_et_look.create_dates(startdate, enddate, 6)

#     cmeta = {
#         "composite_type": "mean",
#         "temporal_interp": False,
#         "spatial_interp": "nearest",
#         "var_name": "ndvi",
#         "var_unit": "-",
#     }

#     dbs = [['/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.06.26.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.07.28.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.07.12.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.05.25.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.08.13.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.06.10.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MOD13/NDVI_MOD13Q1_-_16-daily_2021.08.29.tif'],
#     ['/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.06.18.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.08.21.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.07.04.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.07.20.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.08.05.tif',
#     '/Users/hmcoerver/On My Mac/pyWAPOR/example_data/RAW/MODIS/MYD13/NDVI_MYD13Q1_-_16-daily_2021.06.02.tif']]
    
#     temp_folder = None
#     example_ds = None
#     lean_output = False 
#     diagnostics = None

#     main(cmeta, dbs, epochs_info, lean_output = False)