# -*- coding: utf-8 -*-
"""
Generates input data for `pywapor.se_root`.
"""
from pywapor.collect import downloader
from pywapor.general.logger import log, adjust_logger
from pywapor.general import compositer
import pywapor.general.levels as levels#import pre_et_look_levels, find_example
from pywapor.general import aligner
import datetime

def rename_meteo(ds, *args):
    ds = ds.rename_vars({
                    "t_air": "t_air_i",
                    "u2m": "u2m_i",
                    "v2m": "v2m_i",
                    "qv": "qv_i",
                    "p_air": "p_air_i",
                    "wv": "wv_i",
                    "p_air_0": "p_air_0_i"})
    return ds

def main(folder, latlim, lonlim, timelim, sources, bin_length,
            enhancers = "default", example_source = None, **kwargs):
    """Prepare input data for `se_root`.

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
    example_source : tuple, optional
        Which source to use for spatial alignment, overrides product selected
        through sources, by default None.

    Returns
    -------
    xr.Dataset
        Dataset with data for `pywapor.se_root`.
    """
    t1 = datetime.datetime.now()
    log.info("> PRE_SE_ROOT").add()

    if isinstance(sources, str):
        sources = levels.pre_se_root_levels(sources)

    if isinstance(example_source, type(None)):
        example_source = levels.find_example(sources)
        log.info(f"--> Example dataset is {example_source[0]}.{example_source[1]}.")

    if enhancers == "default":
        enhancers = [rename_meteo]

    bins = compositer.time_bins(timelim, bin_length)
    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])
    ds = aligner.main(dss, sources, example_source, folder, enhancers)

    t2 = datetime.datetime.now()
    log.sub().info(f"< PRE_SE_ROOT ({str(t2 - t1)})")

    return ds

if __name__ == "__main__":

    from pywapor.se_root import main as se_root

    sources = "level_1"

    enhancers = "default"

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 21)]
    bin_length = "DEKAD"
    example_source = None

    sources = levels.pre_se_root_levels(sources)
    example = levels.find_example(sources)

    # print(example)
    ds = main(folder, latlim, lonlim, timelim, sources, bin_length)
    
    ds_out = se_root(ds)