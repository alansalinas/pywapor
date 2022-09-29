# -*- coding: utf-8 -*-
"""
Generates input data for `pywapor.se_root`.
"""
from pywapor.collect import downloader
from pywapor.general import compositer
from pywapor.enhancers.temperature import bt_to_lst
from pywapor.enhancers.other import drop_empty_times
import pywapor.general.levels as levels
from pywapor.general import aligner
import datetime
import os
import numpy as np
import pywapor.general.pre_defaults as defaults
from pywapor.general.logger import log

def rename_meteo(ds, *args):
    """Rename some variables in a dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset whose variables will be renamed.

    Returns
    -------
    xr.Dataset
        Dataset with renamed variables.
    """
    renames = {
                "t_air": "t_air_i",
                "u2m": "u2m_i",
                "v2m": "v2m_i",
                "u": "u_i",
                "qv": "qv_i",
                "p_air": "p_air_i",
                "wv": "wv_i",
                "p_air_0": "p_air_0_i",
                "t_dew": "t_dew_i",
            }
    renames_filtered = {k: v for k,v in renames.items() if k in ds.data_vars}
    ds = ds.rename_vars(renames_filtered)
    return ds

def add_constants(ds, *args):
    """Adds default dimensionless constants to a dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to which the constants will be added.

    Returns
    -------
    xr.Dataset
        Dataset with extra variables.
    """
    ds = ds.assign(defaults.constants_defaults())
    return ds

def main(folder, latlim, lonlim, timelim, sources = "level_1", bin_length = "DEKAD", enhancers = [], **kwargs):
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
    sources : str | dict
        Configuration for each variable and source.
    bin_length : int | "DEKAD"
        Composite length.
    enhancers : list | "default", optional
        Functions to apply to the xr.Dataset before creating the final
        output, by default "default".
    example_source : tuple, optional
        Which source to use for spatial alignment, overrides product selected
        through `sources`, by default None.

    Returns
    -------
    xr.Dataset
        Dataset with input data for `pywapor.se_root`.
    """
    t1 = datetime.datetime.now()
    log.info("> PRE_SE_ROOT").add()

    if not os.path.exists(folder):
        os.makedirs(folder)

    if isinstance(sources, str):
        sources = levels.pre_se_root_levels(sources)

    example_t_vars = [x for x in ["lst", "bt"] if x in sources.keys()]
    example_sources = {k:v for k,v in sources.items() if k in example_t_vars}
    other_sources = {k:v for k,v in sources.items() if k not in example_t_vars}

    general_enhancers = enhancers + [rename_meteo, add_constants, bt_to_lst, drop_empty_times]

    bins = compositer.time_bins(timelim, bin_length)
    adjusted_timelim = [bins[0], bins[-1]]
    buffered_timelim = [adjusted_timelim[0] - np.timedelta64(3, "D"), 
                        adjusted_timelim[1] + np.timedelta64(3, "D")]

    example_dss, example_sources = downloader.collect_sources(folder, example_sources, latlim, lonlim, adjusted_timelim, return_fps = False)

    if len(example_dss) == 0:
        lbl = f"Unable to collect the essential variable(s) (`{'and'.join(example_t_vars)}`) to which the other variables should be aligned."
        log.error("--> " + lbl)
        raise ValueError(lbl)
    
    other_dss, other_sources = downloader.collect_sources(folder, other_sources, latlim, lonlim, buffered_timelim, return_fps = False)
    dss= {**example_dss, **other_dss}

    ds = aligner.main(dss, sources, folder, general_enhancers, example_t_vars = example_t_vars)

    t2 = datetime.datetime.now()
    log.sub().info(f"< PRE_SE_ROOT ({str(t2 - t1)})")

    return ds

if __name__ == "__main__":

    # from pywapor.se_root import main as se_root

    # sources = "level_1"

    enhancers = []

    # folder = r"/Users/hmcoerver/pywapor_notebooks_2"
    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    # timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 7, 11)]
    # bin_length = "DEKAD"
    # example_source = None

    # # _ = adjust_logger(True, folder, "INFO")

    # sources = levels.pre_se_root_levels(sources)

    # sources["ndvi"]["products"] = [
    #     {'source': 'MODIS',
    #         'product_name': 'MOD13Q1.061',
    #         'enhancers': 'default'},
    #     {'source': 'MODIS', 
    #         'product_name': 'MYD13Q1.061', 
    #         'enhancers': 'default'},
    #     {'source': 'PROBAV',
    #         'product_name': 'S5_TOC_100_m_C1',
    #         'enhancers': 'default',
    #         'is_example': True}
    # ]

    # # example = levels.find_example(sources)

    # # print(example)
    # ds = main(folder, latlim, lonlim, timelim, sources, bin_length)
    
    # chunk_size = "20MiB"
    # ds = open_ds(r"/Users/hmcoerver/pywapor_notebooks_2/se_root_in.nc", chunk_size = chunk_size)
    # ds_out = se_root(ds)
    # ds_out = ds_out[["se_root"]]
    # save_ds(ds_out, r"/Users/hmcoerver/pywapor_notebooks_2/test_1.nc", chunk_size = chunk_size)

    # import xarray as xr
    # from dask.diagnostics import ProgressBar

    # ds = xr.open_dataset(r"/Users/hmcoerver/pywapor_notebooks_2/se_root_in.nc", 
    #                         chunks= {"time": -1, "x": 10,"y": 10})
    # ds_out = se_root(ds)
    # ds_out = ds_out[["se_root"]]

    # with ProgressBar():
    #     ds_out.to_netcdf(r"/Users/hmcoerver/pywapor_notebooks_2/test_3.nc")