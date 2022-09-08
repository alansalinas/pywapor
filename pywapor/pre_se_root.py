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
    # TODO remove, this is duplicate from pywapor.pre_et_look
    ds = ds.assign(defaults.constants_defaults())
    return ds

def main(folder, latlim, lonlim, timelim, sources = "level_1", bin_length = "DEKAD",
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

    if isinstance(timelim[0], str):
        timelim[0] = datetime.datetime.strptime(timelim[0], "%Y-%m-%d")
        timelim[1] = datetime.datetime.strptime(timelim[1], "%Y-%m-%d")

    if not os.path.exists(folder):
        os.makedirs(folder)

    if isinstance(sources, str):
        sources = levels.pre_se_root_levels(sources)

    if isinstance(example_source, type(None)):
        example_source = levels.find_example(sources)
        log.info(f"--> Example dataset is {example_source[0]}.{example_source[1]}.")

    example_t_vars = [x for x in ["lst", "bt"] if x in sources.keys()]
    example_sources = {k:v for k,v in sources.items() if k in example_t_vars}
    other_sources = {k:v for k,v in sources.items() if k not in example_t_vars}

    if enhancers == "default":
        enhancers = [rename_meteo,
                    add_constants,
                    bt_to_lst,
                    drop_empty_times,
                    ]

    bins = compositer.time_bins(timelim, bin_length)
    adjusted_timelim = [bins[0], bins[-1]]
    buffered_timelim = [adjusted_timelim[0] - np.timedelta64(3, "D"), 
                        adjusted_timelim[1] + np.timedelta64(3, "D")]

    example_dss = downloader.collect_sources(folder, example_sources, latlim, lonlim, adjusted_timelim)
    other_dss = downloader.collect_sources(folder, other_sources, latlim, lonlim, buffered_timelim)
    dss= {**example_dss, **other_dss}

    ds = aligner.main(dss, sources, example_source, folder, enhancers, example_t_vars = example_t_vars)

    t2 = datetime.datetime.now()
    log.sub().info(f"< PRE_SE_ROOT ({str(t2 - t1)})")

    return ds

# if __name__ == "__main__":

    # from pywapor.se_root import main as se_root

    # sources = "level_1"

    # enhancers = "default"

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