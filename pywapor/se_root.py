import xarray as xr
import os
import numpy as np
import pandas as pd
import pywapor
import pywapor.general.processing_functions as PF
import pywapor.et_look_dev as ETLook_dev
import pywapor.et_look_v2 as ETLook_v2
import tqdm
from pywapor.general.logger import log, adjust_logger
import pywapor.enhancers.dem as dem

def main(level_folder, dss, example_info, et_look_version = "v2"):

    ds_lst, ds_meteo, ds_ndvi = dss
    example_fh, example_ds, example_geoinfo = example_info

    log_write = True
    log_level = "INFO"
    adjust_logger(log_write, level_folder, log_level)
    log.info("> SE_ROOT").add()

    # Version
    if et_look_version == "v2":
        ETLook = ETLook_v2
    elif et_look_version == "dev":
        ETLook = ETLook_dev

    # spatial interpolation
    ds_lst2 = ds_lst.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
    ds_meteo2 = ds_meteo.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)
    ds_ndvi2 = ds_ndvi.interp_like(example_ds, method = "linear", kwargs={"fill_value": "extrapolate"},)

    # temporal interpolation
    ds_meteo3 = ds_meteo2.interp(time = ds_lst2.time, method = "nearest", kwargs={"fill_value": "extrapolate"},) # TODO download more meteo data and switch to linear
    ds_ndvi3 = ds_ndvi2.interp(time = ds_lst.time, method = "linear", kwargs={"fill_value": "extrapolate"},)

    ds_se_root = xr.merge([ds_lst2, ds_meteo3, ds_ndvi3])

    ds_se_root["u_i"] = np.sqrt(ds_se_root["v2m_i"]**2 + ds_se_root["u2m_i"]**2)

    ds_se_root = ds_se_root.rename({"lat": "lat_deg", "lon": "lon_deg"})

    se_root_folder = os.path.join(level_folder, "smc")
    if not os.path.exists(se_root_folder):
        os.mkdir(se_root_folder)

    files = list()

    log.info("--> Calculating se_root.")
    for t in tqdm.tqdm(ds_se_root.time.values): # TODO adjust pywapor.et_look.se_root to remove this forloop

        id = ds_se_root.sel(time = t)
        date = pd.Timestamp(t)

        dt_str = date.strftime("%Y.%m.%d.%H.%M")
        fn = f"SMC_-_-_inst_{dt_str}.tif"
        fh = os.path.join(se_root_folder, fn)

        if not os.path.exists(fh): # TODO also check if fh has correct geoinfo
            se_root_i = pywapor.et_look.se_root(id, None, ETLook, date, verbose = True)[0]["se_root"]
            PF.Save_as_tiff(fh, se_root_i.values, example_geoinfo[0], example_geoinfo[1])

        files.append(fh)

    log.sub().info("< SE_ROOT")

    return [files]