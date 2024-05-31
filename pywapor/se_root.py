# -*- coding: utf-8 -*-
"""
Code to run the SERoot model.
"""
import xarray as xr
import os
import numpy as np
import warnings
import pandas as pd
import datetime
import pywapor.general as g
import pywapor.et_look_v2_v3 as ETLook_v2_v3
from pywapor.general.logger import log, adjust_logger
from pywapor.general.processing_functions import save_ds, open_ds
import copy
import pywapor.pre_se_root as pre_se_root
from pywapor.general import levels

def se_root(folder, latlim, lonlim, timelim, sources = "level_1", bin_length = "DEKAD", se_root_version = None, **kwargs):
    """Runs `pre_se_root` and `se_root`, used internally to generate the `se_root` data in `pre_et_look`.

    Parameters
    ----------
    folder : str
        Path to folder in which to store results.
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

    Returns
    -------
    xr.Dataset
        Dataset with `se_root` variable.
    """
    if isinstance(se_root_version, type(None)) and isinstance(sources, str):
        se_root_version = {True: "v2", False: "v3"}["v3" not in sources]
    elif isinstance(se_root_version, type(None)) and isinstance(sources, dict):
        se_root_version = "v2"
        log.warning("--> Using se_root v2.")

    if isinstance(sources, str):
        sources = levels.pre_se_root_levels(sources)

    ds_in = pre_se_root.main(folder, latlim, lonlim, timelim, sources, bin_length, buffer_timelim = False)
    ds_out = main(ds_in, se_root_version = se_root_version)

    return ds_out

def main(input_data, se_root_version = "v2", export_vars = "default", chunks = {"time": -1, "x": 500, "y": 500}):
    """Run the `se_root` model.

    Parameters
    ----------
    input_data : str | xr.Dataset
        (Path to) dataset generated by `pywapor.pre_se_root`.
    se_root_version : "v2" | "v3", optional
        Which version of the SERoot model to use, by default "v2".
    export_vars : "default" | "all" | list, optional
        Specify which variables to save inside the output file. `"Default"` only 
        stores `se_root`. `"all"` stores all calculated variables. Use a
        list to specify a custom output set, by default "default".
    chunks : dict, optional
        Specify how the calculations are split up. Increase chunk sizes to speed up calculation, 
        decrease to use less RAM, by default {"time": 1, "x": 1000, "y": 1000}.

    Returns
    -------
    xr.Dataset
        Outputs from the model.
    """

    # Inputs
    if isinstance(input_data, str):
        ds = open_ds(input_data, chunks = chunks)
    else:
        ds = input_data.chunk(chunks)
        input_data = ds.encoding["source"]

    _ = adjust_logger(True, os.path.split(input_data)[0], "INFO")

    t1 = datetime.datetime.now()
    log.info("> SE_ROOT").add()

    # Version
    ETLook = ETLook_v2_v3

    log.info(f"--> Running `se_root` ({se_root_version}).").add()

    # Allow skipping of et_look-functions if not all of its required inputs are
    # available.
    g.lazifier.decorate_submods(ETLook, g.lazifier.etlook_decorator)

    ds = g.variables.initiate_ds(ds)

    fp, fn = os.path.split(input_data)

    if not ds["v2m_i"].dtype == object and not ds["u2m_i"].dtype == object:
        ds["u_i"] = np.sqrt(ds["v2m_i"]**2 + ds["u2m_i"]**2)

    doy = [int(pd.Timestamp(x).strftime("%j")) for x in ds["time"].values]
    ds["doy"] = xr.DataArray(doy, coords = ds["time"].coords).chunk("auto")
    dtime = [pd.Timestamp(x).hour + (pd.Timestamp(x).minute / 60) for x in ds["time"].values]
    ds["dtime"] = xr.DataArray(dtime, coords = ds["time"].coords).chunk("auto")
    ds["sc"] = ETLook.solar_radiation.seasonal_correction(ds["doy"])
    ds["decl"] = ETLook.solar_radiation.declination(ds["doy"])
    ds["day_angle"] = ETLook.clear_sky_radiation.day_angle(ds["doy"])

    ds["p_air_i"] = ETLook.meteo.air_pressure_kpa2mbar(ds["p_air_i"])
    ds["p_air_0_i"] = ETLook.meteo.air_pressure_kpa2mbar(ds["p_air_0_i"])

    ds["vc"] = ETLook.leaf.vegetation_cover(ds["ndvi"], nd_min = ds["nd_min"], nd_max = ds["nd_max"], vc_pow = ds["vc_pow"])

    ds["t_air_k_i"] = ETLook.meteo.air_temperature_kelvin_inst(ds["t_air_i"])

    ds["vp_i"] = ETLook.meteo.vapour_pressure_from_specific_humidity_inst(ds["qv_i"], ds["p_air_i"])
    if ds["vp_i"].dtype == object:
        ds["vp_i"] = ETLook.meteo.vapour_pressure_from_dewpoint_inst(ds["t_dew_i"])

    ds["qv_i"] = ETLook.meteo.specific_humidity_from_vapour_pressure(ds["vp_i"], ds["p_air_i"])
    
    ds["ad_moist_i"] = ETLook.meteo.moist_air_density_inst(ds["vp_i"], ds["t_air_k_i"])
    ds["ad_dry_i"] = ETLook.meteo.dry_air_density_inst(ds["p_air_i"], ds["vp_i"], ds["t_air_k_i"])
    ds["ad_i"] = ETLook.meteo.air_density_inst(ds["ad_dry_i"], ds["ad_moist_i"])
    ds["u_b_i_bare"] = ETLook.soil_moisture.wind_speed_blending_height_bare(ds["u_i"], z0m_bare = ds["z0m_bare"], z_obs = ds["z_obs"], z_b = ds["z_b"])
    ds["lon_rad"] = ETLook.solar_radiation.longitude_rad(ds["x"]).chunk("auto")
    ds["lat_rad"] = ETLook.solar_radiation.latitude_rad(ds["y"]).chunk("auto")
    ds["ha"] = ETLook.solar_radiation.hour_angle(ds["sc"], ds["dtime"], ds["lon_rad"])

    ds["ied"] = ETLook.clear_sky_radiation.inverse_earth_sun_distance(ds["day_angle"])
    ds["h0"] = ETLook.clear_sky_radiation.solar_elevation_angle(ds["lat_rad"], ds["decl"], ds["ha"])
    ds["h0ref"] = ETLook.clear_sky_radiation.solar_elevation_angle_refracted(ds["h0"])
    ds["m"] = ETLook.clear_sky_radiation.relative_optical_airmass(ds["p_air_i"], ds["p_air_0_i"], ds["h0ref"])
    ds["rotm"] = ETLook.clear_sky_radiation.rayleigh_optical_thickness(ds["m"])
    ds["Tl2"] = ETLook.clear_sky_radiation.linke_turbidity(ds["wv_i"], ds["aod550_i"], ds["p_air_i"], ds["p_air_0_i"])
    ds["G0"] = ETLook.clear_sky_radiation.extraterrestrial_irradiance_normal(ds["IO"], ds["ied"])
    ds["B0c"] = ETLook.clear_sky_radiation.beam_irradiance_normal_clear(ds["G0"], ds["Tl2"], ds["m"], ds["rotm"], ds["h0"])
    ds["Bhc"] = ETLook.clear_sky_radiation.beam_irradiance_horizontal_clear(ds["B0c"], ds["h0"])
    ds["Dhc"] = ETLook.clear_sky_radiation.diffuse_irradiance_horizontal_clear(ds["G0"], ds["Tl2"], ds["h0"])

    ds["ra_hor_clear_i"] = ETLook.clear_sky_radiation.ra_clear_horizontal(ds["Bhc"], ds["Dhc"])
    ds["emiss_atm_i"] = ETLook.soil_moisture.atmospheric_emissivity_inst(ds["vp_i"], ds["t_air_k_i"])

    ds["rn_bare"] = ETLook.soil_moisture.net_radiation_bare(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["lst"], ds["r0_bare"])
    ds["rn_full"] = ETLook.soil_moisture.net_radiation_full(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["lst"], ds["r0_full"])

    ds["h_bare"] = ETLook.soil_moisture.sensible_heat_flux_bare(ds["rn_bare"], fraction_h_bare = ds["fraction_h_bare"])
    ds["h_full"] = ETLook.soil_moisture.sensible_heat_flux_full(ds["rn_full"], fraction_h_full = ds["fraction_h_full"])
    ds["u_b_i_full"] = ETLook.soil_moisture.wind_speed_blending_height_full_inst(ds["u_i"], z0m_full = ds["z0m_full"], z_obs = ds["z_obs"], z_b = ds["z_b"])

    ds["u_star_i_bare"] = ETLook.soil_moisture.friction_velocity_bare_inst(ds["u_b_i_bare"], z0m_bare = ds["z0m_bare"], disp_bare = ds["disp_bare"], z_b = ds["z_b"])
    ds["u_star_i_full"] = ETLook.soil_moisture.friction_velocity_full_inst(ds["u_b_i_full"], z0m_full = ds["z0m_full"], disp_full = ds["disp_full"], z_b = ds["z_b"])
    ds["L_bare"] = ETLook.soil_moisture.monin_obukhov_length_bare(ds["h_bare"], ds["ad_i"], ds["u_star_i_bare"], ds["t_air_k_i"])
    ds["L_full"] = ETLook.soil_moisture.monin_obukhov_length_full(ds["h_full"], ds["ad_i"], ds["u_star_i_full"], ds["t_air_k_i"])

    ds["u_i_soil"] = ETLook.soil_moisture.wind_speed_soil_inst(ds["u_i"], ds["L_bare"], z_obs = ds["z_obs"])
    
    ds["ras"] = ETLook.soil_moisture.aerodynamical_resistance_forced_convection_soil(ds["u_i_soil"])
    ds["raa"] = ETLook.soil_moisture.aerodynamical_resistance_forced_convection_bare(ds["u_i"], ds["L_bare"], z0m_bare = ds["z0m_bare"], disp_bare = ds["disp_bare"], z_obs = ds["z_obs"])
    ds["rac"] = ETLook.soil_moisture.aerodynamical_resistance_forced_convection_full(ds["u_i"], ds["L_full"], z0m_full = ds["z0m_full"], disp_full = ds["disp_full"], z_obs = ds["z_obs"])
    
    if se_root_version == "v3":

        ds["rah_bare_free"] = ETLook.soil_moisture.aerodynamical_resistance_free_convection_bare(ds["h_bare"], ds["t_air_k_i"], ds["ad_i"], z0m_bare = ds["z0m_bare"])
        ds["rah_full_free"] = ETLook.soil_moisture.aerodynamical_resistance_free_convection_full(ds["h_full"], ds["t_air_k_i"], ds["ad_i"], z0m_full = ds["z0m_full"])
        
        ds["raa"] = ETLook.soil_moisture.aerodynamical_resistance_bare(ds["raa"], ds["ras"], rah_bare_free = ds["rah_bare_free"])
        ds["rac"] = ETLook.soil_moisture.aerodynamical_resistance_full(ds["rac"], rah_full_free = ds["rah_full_free"])

    ds["t_max_bare"] = ETLook.soil_moisture.maximum_temperature_bare(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["raa"], ds["ras"], ds["r0_bare"])
    ds["t_max_full"] = ETLook.soil_moisture.maximum_temperature_full(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["rac"], ds["r0_full"])

    ds["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst_new(ds["t_air_i"], ds["qv_i"], ds["p_air_i"])
    ds["lst_max"] = ETLook.soil_moisture.maximum_temperature(ds["t_max_bare"], ds["t_max_full"], ds["vc"])

    ds["t_wet_k_i"] = ETLook.meteo.wet_bulb_temperature_kelvin_inst(ds["t_wet_i"])
    ds["lst_min"] = ETLook.soil_moisture.minimum_temperature(ds["t_wet_k_i"], ds["t_air_k_i"], ds["vc"])

    ds["se_root"] = ETLook.soil_moisture.soil_moisture_from_maximum_temperature(ds["lst_max"], ds["lst"], ds["lst_min"])

    if export_vars == "all":
        ...
    elif export_vars == "default":
        keep_vars = [
            'se_root'
            ]
        keep_vars = [x for x in keep_vars if x in ds.variables]
        ds = ds[keep_vars]
    elif isinstance(export_vars, list):
        keep_vars = copy.copy(export_vars)
        ds = ds[keep_vars]
    else:
        raise ValueError(f"Please provide a valid `export_vars` ('all', 'default' or a list).")

    ds = g.variables.fill_attrs(ds)

    log.sub()
    if len(ds.data_vars) == 0:
        log.info("--> No data to export, try adjusting `export_vars`.")
        ds = None
    else:
        fn = fn.replace("in", "out")
        fp_out = os.path.join(fp, fn)
        while os.path.isfile(fp_out):
            fp_out = fp_out.replace(".nc", "_.nc")

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in power")
            warnings.filterwarnings("ignore", message="invalid value encountered in log")
            ds = save_ds(ds, fp_out, encoding = "initiate", chunks = chunks, label = f"Saving output to `{os.path.split(fp_out)[-1]}`.")

    t2 = datetime.datetime.now()
    log.sub().info(f"< SE_ROOT ({str(t2 - t1)})")

    return ds

def test_ds(ds, var):
    from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler, ProgressBar
    from dask.diagnostics import visualize

    with Profiler() as prof, ResourceProfiler(dt=0.25) as rprof, CacheProfiler() as cprof, ProgressBar():
        out = ds[var].compute()

    visualize([prof, rprof, cprof])
    return out

if __name__ == "__main__":

    se_root_version = "v2"
    export_vars = "default"
    chunks = {"time": 1, "x": 1000, "y": 1000}
    # input_data = r"/Users/hmcoerver/Local/20220325_20220415_test_data/se_root_in.nc"

    # out = main(input_data, se_root_version = se_root_version, export_vars = export_vars)
