from pywapor.general.processing_functions import save_ds, open_ds
from pywapor.general.reproject import reproject
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.logger import log
import os
import xarray as xr

def main(dss, sources, example_source, folder, enhancers, example_t_var = "lst"):
    """Aligns the datetimes in de `dss` xr.Datasets with the datetimes of the 
    dataset with variable `example_t_var`.

    Parameters
    ----------
    dss : dict
        Keys are tuples of (`source`, `product_name`), values are xr.Dataset's 
        which will be aligned along the time dimensions.
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
    final_path = os.path.join(folder, "se_root_in.nc")
    if os.path.isfile(final_path):
        return open_ds(final_path, "all")

    example_ds = open_ds(dss[example_source], "all")

    dss2 = list()

    # Move example_t_var to first position.
    variables = list(sources.keys())
    variables.insert(0, variables.pop(variables.index(example_t_var)))

    for var in variables:

        config = sources[var]
        spatial_interp = config["spatial_interp"]
        temporal_interp = config["temporal_interp"]
        srcs = [(x["source"], x["product_name"]) for x in config["products"]]

        dst_path = os.path.join(folder, f"{var}_i.nc")
        if os.path.isfile(dst_path):
            ds = open_ds(dst_path, "all")
            if var == example_t_var:
                example_time = ds["time"]
            dss2.append(ds)
            continue

        # Align pixels.
        dss1 = [reproject(open_ds(dss[src])[[var]], example_ds, dst_path.replace(".nc", f"_{i}.nc"), spatial_interp = spatial_interp) for i, src in enumerate(srcs)]

        # Combine different source_products (along time dimension).
        ds = xr.combine_nested(dss1, concat_dim = "time").chunk({"time": -1}).sortby("time").squeeze()

        if var == example_t_var:
            example_time = ds["time"]
        elif "time" in ds[var].dims:
            log.info(f"--> Aligning times in `{var}` with `{example_t_var}` ({temporal_interp}).")
            ds = ds.interp_like(example_time, method = temporal_interp)
            ds = ds.ffill("time").bfill("time") # NOTE extrapolate.

        # Save output
        dss2.append(save_ds(ds, dst_path, decode_coords = "all"))

        for nc in dss1:
            os.remove(nc.encoding["source"])

    ds = xr.merge(dss2)

    # Apply product specific functions.
    for func in enhancers:
        ds, label = apply_enhancer(ds, var, func)
        log.info(label)

    ds = ds.drop_vars([x for x in ds.coords if (x not in ds.dims) and (x != "spatial_ref")])
    ds = save_ds(ds, final_path, decode_coords = "all")

    for nc in dss2:
        os.remove(nc.encoding["source"])

    return ds

if __name__ == "__main__":

    example_t_var = "lst"

#     import datetime
#     from pywapor.general import compositer
#     from pywapor.collect import downloader

#     sources = {
#         "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
#         "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
#         "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
#         "u2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
#         "v2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
#         "qv":           [("MERRA2", "M2I1NXASM.5.12.4")],
#         "wv":           [("MERRA2", "M2I1NXASM.5.12.4")],
#         "p_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
#         "p_air_0":      [("MERRA2", "M2I1NXASM.5.12.4")],
#     }

#     folder = r"/Users/hmcoerver/Downloads/pywapor_test"
#     latlim = [28.9, 29.7]
#     lonlim = [30.2, 31.2]
#     timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
#     example_source = ("MODIS", "MOD13Q1.061")
#     bin_length = 4

#     bins = compositer.time_bins(timelim, bin_length)
#     dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

#     main(dss, sources, example_source, folder, example_t_var = "lst")