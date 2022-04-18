import pywapor
from pywapor.general.processing_functions import save_ds, open_ds
from pywapor.general.reproject import reproject
import os
import xarray as xr

def main(dss, sources, example_source, folder, example_t_var = "lst"):

    example_ds = open_ds(dss[example_source], "all")

    cdefaults = pywapor.general.pre_defaults.composite_defaults()

    dss2 = list()

    # Move example_t_var to first position.
    variables = list(sources.keys())
    variables.insert(0, variables.pop(variables.index(example_t_var)))

    for var in variables:

        srcs = sources[var]

        spatial_interp = cdefaults[var]["spatial_interp"]
        temporal_interp = cdefaults[var]["temporal_interp"]

        dst_path = os.path.join(folder, f"{var}.nc")

        # Align pixels.
        dss1 = [reproject(open_ds(dss[src])[[var]], example_ds, dst_path.replace(".nc", f"_{i}.nc"), spatial_interp = spatial_interp) for i, src in enumerate(srcs)]

        # Combine different source_products (along time dimension).
        ds = xr.combine_nested(dss1, concat_dim = "time").sortby("time").squeeze()

        if (var == example_t_var) and ("time" in ds[var].dims):
            example_time = ds["time"]

        elif ("time" in ds[var].dims):
            # ds = ds.chunk({"time": -1})
            ds = ds.interp_like(example_time, method = temporal_interp)

        # Save output
        dss2.append(save_ds(ds, dst_path, decode_coords = "all"))

        for nc in dss1:
            os.remove(nc.encoding["source"])

    final_path = os.path.join(folder, "se_root_in.nc")

    ds = xr.merge(dss2)
    ds = ds.drop_vars([x for x in ds.coords if (x not in ds.dims) and (x != "spatial_ref")])
    ds = save_ds(ds, final_path, decode_coords = "all")

    for nc in dss2:
        os.remove(nc.encoding["source"])

    return ds

if __name__ == "__main__":

    import datetime
    from pywapor.general import compositer
    from pywapor.collect import downloader

    sources = {

        "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
        "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
        "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        "u2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        "v2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        "qv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        "wv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        "p_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        "p_air_0":      [("MERRA2", "M2I1NXASM.5.12.4")],

    }

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
    example_source = ("MODIS", "MOD13Q1.061")
    bin_length = 4

    bins = compositer.time_bins(timelim, bin_length)
    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

    main(dss, sources, example_source, folder, example_t_var = "lst")