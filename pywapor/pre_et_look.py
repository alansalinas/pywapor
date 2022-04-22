"""
Generates input data for `pywapor.et_look`.
"""
from pywapor.collect import downloader
from pywapor.general.logger import log, adjust_logger
from pywapor.general import compositer
import pywapor.general.levels as levels
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from functools import partial
import pywapor.general.pre_defaults as defaults
from pywapor.general.variables import fill_attrs
from pywapor.enhancers.temperature import lapse_rate as _lapse_rate

def rename_vars(ds, *args):
    varis = ["p", "ra", "t_air", "t_air_min", "t_air_max", 
            "u2m", "v2m", "qv", "p_air", "p_air_0", "wv"]
    present_vars = [x for x in varis if x in ds.variables]
    ds = ds.rename({k: k + "_24" for k in present_vars})
    return ds

def lapse_rate(ds, *args):
    present_vars = [x for x in ds.variables if "t_air" in x]
    for var in present_vars:
        ds = _lapse_rate(ds, var)
    return ds

def calc_doys(ds, *args, bins = None):
    bin_doys = [int(pd.Timestamp(x).strftime("%j")) for x in bins]
    doy = np.mean([bin_doys[:-1], bin_doys[1:]], axis=0, dtype = int)
    ds["doy"] = xr.DataArray(doy, coords = ds["time_bins"].coords).chunk("auto")
    return ds

def add_constants(ds, *args):
    ds = ds.assign(defaults.constants_defaults())
    return ds

def main(folder, latlim, lonlim, timelim, sources, bin_length, 
            enhancers = "default", diagnostics = None, example_source = None):
    """Prepare input data for `et_look`.

    Parameters
    ----------
    folder : str
        Path to folder in which to store (intermediate) data.
    latlim : list
        Latitude limits of area of interest.
    lonlim : list
        Longitude limits of area of interest.
    timelim : list
        Period for which to prepare data.
    sources : dict
        Configuration for each variable and source.
    bin_length : int | "DEKAD"
        Composite length.
    enhancers : list | "default", optional
        Functions to apply to the xr.Dataset before creating the final
        output, by default "default".
    diagnostics : dict, optional
        Dictionary with coordinates and point-labels for which graphs can be 
        created.
    example_source : tuple, optional
        Which source to use for spatial alignment, overrides product selected
        through sources, by default None.

    Returns
    -------
    xr.Dataset
        Dataset with data for `pywapor.et_look`.
    """

    _ = adjust_logger(True, folder, "INFO")

    t1 = datetime.datetime.now()
    log.info("> PRE_ET_LOOK").add()

    if isinstance(sources, str):
        sources = levels.pre_et_look_levels(sources)

    if isinstance(example_source, type(None)):
        example_source = levels.find_example(sources)
        log.info(f"--> Example dataset is {example_source[0]}.{example_source[1]}.")

    bins = compositer.time_bins(timelim, bin_length)

    if enhancers == "default":
        enhancers = [
                    rename_vars, 
                    fill_attrs, 
                    lapse_rate, 
                    partial(calc_doys, bins = bins),
                    add_constants
                    ]

    # sources.pop("se_root")
    sources = {k:v for k, v in sources.items() if k in ["ndvi"]}

    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

    if diagnostics:
        t_1 = datetime.datetime.now()
        log.info("> DIAGNOSTICS").add()
        ds = compositer.main(dss, sources, example_source, bins, folder, enhancers, diagnostics = diagnostics)
        t_2 = datetime.datetime.now()
        log.sub().info(f"< DIAGNOSTICS ({str(t_2 - t_1)})")
    else:
        ds = compositer.main(dss, sources, example_source, bins, folder, enhancers, diagnostics = None)

    t2 = datetime.datetime.now()
    log.sub().info(f"< PRE_ET_LOOK ({str(t2 - t1)})")

    return ds

if __name__ == "__main__":

    from pywapor.et_look import main as et_look

    diagnostics = { # label          # lat      # lon
                    "water":	    (29.44977,	30.58215),
                    "desert":	    (29.12343,	30.51222),
                    "agriculture":	(29.32301,	30.77599),
                    "urban":	    (29.30962,	30.84109),
                    }

    # diagnostics = None
    example_source = None
    enhancers = "default"

    sources = "level_1"
    sources = levels.pre_et_look_levels(sources)
    sources = {k:v for k, v in sources.items() if k in ["ndvi"]}

    sources["ndvi"]["temporal_interp"] = False
    sources["ndvi"]["composite_type"] = "max"

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 8, 5)]
    # bin_length = "DEKAD"
    bin_length = 15

    ds = main(folder, latlim, lonlim, timelim, sources, 
                bin_length, diagnostics = diagnostics)

    # from pywapor.general.processing_functions import open_ds
    # input_data = open_ds(r"/Users/hmcoerver/Downloads/pywapor_test/et_look_in.nc")
    # ds_out = et_look(ds, export_vars = "default")