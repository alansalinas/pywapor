# -*- coding: utf-8 -*-
"""
Generates input data for `pywapor.et_look`.
"""
from pywapor.collect import downloader
from pywapor.general.logger import log, adjust_logger
from pywapor.general import compositer
import pywapor.general.levels as levels
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

    log.info("> PRE_ET_LOOK").add()

    if isinstance(sources, str):
        sources = levels.pre_et_look_levels(sources)

    if isinstance(example_source, type(None)):
        example_source = levels.find_example(sources)
        log.info(f"--> Example dataset is {example_source[0]}.{example_source[1]}.")

    if enhancers == "default":
        enhancers = [rename_vars, fill_attrs, lapse_rate]

    sources = {k:v for k, v in sources.items() if k in ["ndvi", "t_air_min", "t_air", "z"]}

    bins = compositer.time_bins(timelim, bin_length)
    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])

    if diagnostics:
        log.info("> DIAGNOSTICS").add()
        _ = compositer.main(dss, sources, example_source, bins, folder, enhancers, diagnostics = diagnostics)
        log.sub().info("< DIAGNOSTICS")

    ds = compositer.main(dss, sources, example_source, bins, folder, enhancers, diagnostics = None)

    log.sub().info("< PRE_ET_LOOK")

    return ds

if __name__ == "__main__":

    import datetime

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

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 21)]
    bin_length = "DEKAD"

    ds = main(folder, latlim, lonlim, timelim, sources, 
                bin_length, diagnostics = diagnostics)

    # from pywapor.general.processing_functions import open_ds
    # ds = open_ds(r"/Users/hmcoerver/Downloads/pywapor_test/et_look_in.nc")
    # from pywapor.et_look import main as et_look
    # ds_out = et_look(ds)