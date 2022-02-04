import os
import numpy as np
import pywapor.et_look_dev as ETLook_dev
import pywapor.et_look_v2 as ETLook_v2
import pywapor.general as g
import pywapor.general.processing_functions as PF
import xarray as xr
import pandas as pd
from pywapor.general.compositer import calculate_ds
import copy

def main(input_data, et_look_version = "v2", export_vars = "default", export_to_tif = False):

    # Version
    if et_look_version == "v2":
        ETLook = ETLook_v2
        print("--> Running ETLook_v2")
    elif et_look_version == "dev":
        ETLook = ETLook_dev
        print("--> Running ETLOOK_dev")
    c = ETLook.constants

    # Inputs
    if isinstance(input_data, str):
        ds = xr.open_dataset(input_data)
    else:
        ds = input_data
        input_data = ds.encoding["source"]

    # # spatial_or_constant = ["lw_offset", "lw_slope", "r0_bare", "r0_full", 
    # #                         "rn_offset", "rn_slope", "t_opt", "vpd_slope"]
    # for param in spatial_or_constant:
    #     if isinstance(id[param], type(None)):
    #         id[param] = np.ones_like(id["ndvi"]) * getattr(c, param)

    ds["u_24"] = np.sqrt(ds["v2m_24"]**2 + ds["u2m_24"]**2)

    # Constants
    doy_epoch_start = [int(pd.Timestamp(x).strftime("%j")) for x in ds["epoch_starts"].values]
    doy_epoch_end = [int(pd.Timestamp(x).strftime("%j")) for x in ds["epoch_ends"].values]
    doy = [int((x+y)/2) for x, y in zip(doy_epoch_start, doy_epoch_end)]
    ds["doy"] = xr.DataArray(doy, coords = ds["epoch_starts"].coords)

    ds["sc"] = ETLook.solar_radiation.seasonal_correction(ds["doy"])
    ds["decl"] = ETLook.solar_radiation.declination(ds["doy"])
    ds["iesd"] = ETLook.solar_radiation.inverse_earth_sun_distance(ds["doy"])

    ######################## MODEL ETLOOK ####################################

    # **effective_leaf_area_index*********************************************
    ds["vc"] = ETLook.leaf.vegetation_cover(ds["ndvi"], c.nd_min, c.nd_max, c.vc_pow)
    ds["lai"] = ETLook.leaf.leaf_area_index(ds["vc"], c.vc_min, c.vc_max, c.lai_pow)
    ds["lai_eff"] = ETLook.leaf.effective_leaf_area_index(ds["lai"])
    
    # *******TRANSPIRATION COMPONENT******************************************
    
    # **soil fraction*********************************************************
    ds["sf_soil"] = ETLook.radiation.soil_fraction(ds["lai"])

    # **atmospheric canopy resistance***********************************************
    ds["lat_rad"] = ETLook.solar_radiation.latitude_rad(ds["lat_deg"])
    ds["ws"] = ETLook.solar_radiation.sunset_hour_angle(ds["lat_rad"], ds["decl"])
    
    if "ra_24" not in list(ds.variables):
        ds["slope_rad"] = ETLook.solar_radiation.slope_rad(ds["slope"]) # TODO adjust pre_et_look to make this called 'slope_deg' again.
        ds["aspect_rad"] = ETLook.solar_radiation.aspect_rad(ds["aspect"]) # TODO adjust pre_et_look to make this called 'aspect_deg' again.
        ds["ra_24_toa"] = ETLook.solar_radiation.daily_solar_radiation_toa_new(ds["sc"], ds["decl"], ds["iesd"], ds["lat_rad"], ds["doy"], ds["slope_rad"], ds["aspect_rad"])
        ds["ra_24_toa_flat"] = ETLook.solar_radiation.daily_solar_radiation_toa_flat(ds["decl"], ds["iesd"], ds["lat_rad"], ds["ws"])
        ds["diffusion_index"] = ETLook.solar_radiation.diffusion_index(ds["trans_24"], c.diffusion_slope, c.diffusion_intercept)
        ds["ra_24"] = ETLook.solar_radiation.daily_total_solar_radiation(ds["ra_24_toa"], ds["ra_24_toa_flat"], ds["diffusion_index"], ds["trans_24"])
    else:
        ds["ra_24_toa_flat"] = ETLook.solar_radiation.daily_solar_radiation_toa_flat(ds["decl"], ds["iesd"], ds["lat_rad"], ds["ws"])
        ds["trans_24"] = ds["ra_24"] / ds["ra_24_toa_flat"]

    ds["stress_rad"] = ETLook.stress.stress_radiation(ds["ra_24"])
    ds["p_air_0_24"] = ETLook.meteo.air_pressure_kpa2mbar(ds["p_air_0_24"])
    ds["p_air_24"] = ETLook.meteo.air_pressure_daily(ds["z"], ds["p_air_0_24"])
    ds["vp_24"] = ETLook.meteo.vapour_pressure_from_specific_humidity_daily(ds["qv_24"], ds["p_air_24"])

    if "t_air_24" not in list(ds.variables):
        ds["t_air_24"] = (ds["t_air_min_24"] + ds["t_air_max_24"]) / 2

    if et_look_version == "v2":
        ds["svp_24"] = ETLook.meteo.saturated_vapour_pressure(ds["t_air_24"])
    elif et_look_version == "dev":
        ds["svp_24"] = ETLook.meteo.saturated_vapour_pressure_average(
                    ETLook.meteo.saturated_vapour_pressure_maximum(ds["t_air_max_24"]),
                    ETLook.meteo.saturated_vapour_pressure_minimum(ds["t_air_min_24"]))

    ds["vpd_24"] = ETLook.meteo.vapour_pressure_deficit_daily(ds["svp_24"], ds["vp_24"])
    ds["stress_vpd"] = ETLook.stress.stress_vpd(ds["vpd_24"], ds["vpd_slope"])
    ds["stress_temp"] = ETLook.stress.stress_temperature(ds["t_air_24"], ds["t_opt"], c.t_min, c.t_max)

    if et_look_version == "dev":
        if isinstance(ds["land_mask"], xr.DataArray):
            ds["rs_min"] = xr.where(ds["land_mask"] == 3, 400, 100)
        else:
            ds["rs_min"] = np.where(ds["land_mask"] == 3, 400, 100)

    ds["r_canopy_0"] = ETLook.resistance.atmospheric_canopy_resistance(ds["lai_eff"], ds["stress_rad"], ds["stress_vpd"], ds["stress_temp"], ds["rs_min"], c.rcan_max)

    # **net radiation canopy******************************************************
    ds["t_air_k_24"] = ETLook.meteo.air_temperature_kelvin_daily(ds["t_air_24"])
    ds["l_net"] = ETLook.radiation.longwave_radiation_fao(ds["t_air_k_24"], ds["vp_24"], ds["trans_24"], c.vp_slope, c.vp_offset, ds["lw_slope"], ds["lw_offset"])
    ds["int_mm"] = ETLook.evapotranspiration.interception_mm(ds["p_24"], ds["vc"], ds["lai"], c.int_max)
    ds["lh_24"] = ETLook.meteo.latent_heat_daily(ds["t_air_24"])
    ds["int_wm2"] = ETLook.radiation.interception_wm2(ds["int_mm"], ds["lh_24"])
    ds["rn_24"] = ETLook.radiation.net_radiation(ds["r0"], ds["ra_24"], ds["l_net"], ds["int_wm2"])
    ds["rn_24_canopy"] = ETLook.radiation.net_radiation_canopy(ds["rn_24"], ds["sf_soil"])

    # find water region using ndvi
    if et_look_version == "dev":
        if ds.pixel_size < 250:
            if isinstance(ds["land_mask"], xr.DataArray):
                ds["land_mask"] = xr.where(ds["land_mask"] == 2, 1, ds["land_mask"])
                ds["land_mask"] = xr.where(ds["ndvi"] < 0, 2, ds["land_mask"])
            else:
                ds["land_mask"] = np.where(ds["land_mask"] == 2, 1, ds["land_mask"])
                ds["land_mask"] = np.where(ds["ndvi"] < 0, 2, ds["land_mask"])

    ds["stress_moist"] = ETLook.stress.stress_moisture(ds["se_root"], c.tenacity)
    ds["r_canopy"] = ETLook.resistance.canopy_resistance(ds["r_canopy_0"], ds["stress_moist"], c.rcan_max)

    # **initial canopy aerodynamic resistance***********************************************************
    
    if "z_oro" not in list(ds.variables): # TODO check if this works
        ds["slope_rad"] = ETLook.solar_radiation.slope_rad(id["slope_deg"])
        ds["z_oro"] = ETLook.roughness.orographic_roughness(ds["slope_rad"], ds.pixel_size)

    ds["z_obst"] = ETLook.roughness.obstacle_height(ds["ndvi"], ds["z_obst_max"], c.ndvi_obs_min, c.ndvi_obs_max, c.obs_fr)
    ds["z0m"] = ETLook.roughness.roughness_length(ds["lai"], ds["z_oro"], ds["z_obst"], ds["z_obst_max"], ds["land_mask"])
    ds["ra_canopy_init"] = ETLook.neutral.initial_canopy_aerodynamic_resistance(ds["u_24"], ds["z0m"], c.z_obs)

    # **windspeed blending height daily***********************************************************
    ds["u_b_24"] = ETLook.meteo.wind_speed_blending_height_daily(ds["u_24"], c.z_obs, c.z_b)

    # **ETLook.neutral.initial_daily_transpiration***********************************************************
    ds["ad_dry_24"] = ETLook.meteo.dry_air_density_daily(ds["p_air_24"], ds["vp_24"], ds["t_air_k_24"])
    ds["ad_moist_24"] = ETLook.meteo.moist_air_density_daily(ds["vp_24"], ds["t_air_k_24"])
    ds["ad_24"] = ETLook.meteo.air_density_daily(ds["ad_dry_24"], ds["ad_moist_24"])
    ds["psy_24"] = ETLook.meteo.psychrometric_constant_daily(ds["p_air_24"], ds["lh_24"])
    ds["ssvp_24"] = ETLook.meteo.slope_saturated_vapour_pressure_daily(ds["t_air_24"])
    ds["t_24_init"] = ETLook.neutral.initial_daily_transpiration(ds["rn_24_canopy"], ds["ssvp_24"], ds["ad_24"], ds["vpd_24"], ds["psy_24"], ds["r_canopy"], ds["ra_canopy_init"])

    # **ETLook.unstable.initial_sensible_heat_flux_canopy_daily***********************************************************
    ds["h_canopy_24_init"] = ETLook.unstable.initial_sensible_heat_flux_canopy_daily(ds["rn_24_canopy"], ds["t_24_init"])

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    ds["disp"] = ETLook.roughness.displacement_height(ds["lai"], ds["z_obst"], ds["land_mask"], c.c1)
    ds["u_star_24_init"] = ETLook.unstable.initial_friction_velocity_daily(ds["u_b_24"], ds["z0m"], ds["disp"], c.z_b)

    # **ETLook.unstable.transpiration***********************************************************
    ds["t_24"] = ETLook.unstable.transpiration(ds["rn_24_canopy"], ds["ssvp_24"], ds["ad_24"], ds["vpd_24"], ds["psy_24"], ds["r_canopy"], ds["h_canopy_24_init"], ds["t_air_k_24"], ds["u_star_24_init"], ds["z0m"], ds["disp"], ds["u_b_24"], c.z_obs, c.z_b, c.iter_h)
    ds["t_24_mm"] = ETLook.unstable.transpiration_mm(ds["t_24"], ds["lh_24"])

    if et_look_version == "dev":
        ds["tpot_24"] = ETLook.unstable.transpiration(ds["rn_24_canopy"], ds["ssvp_24"], ds["ad_24"], ds["vpd_24"], ds["psy_24"], ds["r_canopy"] * ds["stress_moist"], ds["h_canopy_24_init"], ds["t_air_k_24"], ds["u_star_24_init"], ds["z0m"], ds["disp"], ds["u_b_24"], c.z_obs, c.z_b, c.iter_h)
        ds["tpot_24_mm"] = ETLook.unstable.transpiration_mm(ds["tpot_24"], ds["lh_24"])

    #*******EVAPORATION COMPONENT****************************************************************

    # **ETLook.radiation.net_radiation_soil***********************************************************
    ds["rn_24_soil"] = ETLook.radiation.net_radiation_soil(ds["rn_24"], ds["sf_soil"])

    # **ETLook.resistance.soil_resistance***********************************************************
    ds["r_soil"] = ETLook.resistance.soil_resistance(ds["se_root"], ds["land_mask"], c.r_soil_pow, c.r_soil_min)

    # **ETLook.resistance.soil_resistance***********************************************************
    ds["ra_soil_init"] = ETLook.neutral.initial_soil_aerodynamic_resistance(ds["u_24"], c.z_obs)

    # **ETLook.unstable.initial_friction_velocity_soil_daily***********************************************************
    ds["u_star_24_soil_init"] = ETLook.unstable.initial_friction_velocity_soil_daily(ds["u_b_24"], ds["disp"], c.z_b)

    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************
    if et_look_version == "dev":
        ds["stc"] = ETLook.radiation.soil_thermal_conductivity(c.se_top)
        ds["vhc"] = ETLook.radiation.volumetric_heat_capacity(c.se_top, c.porosity)
    elif et_look_version == "v2":
        ds["stc"] = ETLook.radiation.soil_thermal_conductivity(ds["se_root"])
        ds["vhc"] = ETLook.radiation.volumetric_heat_capacity(ds["se_root"], c.porosity)
       
    ds["dd"] = ETLook.radiation.damping_depth(ds["stc"], ds["vhc"])
    ds["g0_bs"] = ETLook.radiation.bare_soil_heat_flux(ds["doy"], ds["dd"], ds["stc"], ds["t_amp_year"], ds["lat_rad"])
    ds["g0_24"] = ETLook.radiation.soil_heat_flux(ds["g0_bs"], ds["sf_soil"], ds["land_mask"], ds["rn_24_soil"], ds["trans_24"], ds["ra_24"], ds["l_net"], ds["rn_slope"], ds["rn_offset"])
    ds["e_24_init"] = ETLook.neutral.initial_daily_evaporation(ds["rn_24_soil"], ds["g0_24"], ds["ssvp_24"], ds["ad_24"], ds["vpd_24"], ds["psy_24"], ds["r_soil"], ds["ra_soil_init"])
    ds["h_soil_24_init"] = ETLook.unstable.initial_sensible_heat_flux_soil_daily(ds["rn_24_soil"], ds["e_24_init"], ds["g0_24"])

    # **ETLook.unstable.evaporation***********************************************************
    ds["e_24"] = ETLook.unstable.evaporation(ds["rn_24_soil"], ds["g0_24"], ds["ssvp_24"], ds["ad_24"], ds["vpd_24"], ds["psy_24"], ds["r_soil"], ds["h_soil_24_init"], ds["t_air_k_24"], ds["u_star_24_soil_init"], ds["disp"], ds["u_b_24"], c.z_b, c.z_obs, c.iter_h)
    ds["e_24_mm"] = ETLook.unstable.evaporation_mm(ds["e_24"], ds["lh_24"])
    ds["et_24_mm"] = ETLook.evapotranspiration.et_actual_mm(ds["e_24_mm"], ds["t_24_mm"])
    
    if et_look_version == "dev":
        ds["e_24_mm"] = np.clip(ds["e_24_mm"], 0, np.inf)
        ds["et_24_mm"] = np.clip(ds["et_24_mm"], 0, np.inf)

    # **ETLook.unstable.evaporation***********************************************************
    ds["rn_24_grass"] = ETLook.radiation.net_radiation_grass(ds["ra_24"], ds["l_net"], c.r0_grass)
    ds["et_ref_24"] = ETLook.evapotranspiration.et_reference(ds["rn_24_grass"], ds["ad_24"], ds["psy_24"], ds["vpd_24"], ds["ssvp_24"], ds["u_24"])
    ds["et_ref_24_mm"] = ETLook.evapotranspiration.et_reference_mm(ds["et_ref_24"], ds["lh_24"])

    if et_look_version == "v2":
        ds["et_ref_24_mm"] = np.clip(ds["et_ref_24_mm"], 0, np.inf)

    if et_look_version == "dev":
        eps_a = ETLook.stress.epsilon_autotrophic_respiration()     
        ds["lue"] = ETLook.biomass.lue(ds["lue_max"], ds["stress_temp"], ds["stress_moist"], eps_a)
        ds["fpar"] = ETLook.leaf.fpar(ds["vc"], ds["ndvi"])
        ds["apar"] = ETLook.leaf.apar(ds["ra_24"], ds["fpar"])       
        ds["biomass_prod"] = ETLook.biomass.biomass(ds["apar"], ds["lue"])     

    if export_vars == "all":
        ...
    elif export_vars == "default":
        keep_vars = ['int_mm',
                    't_24_mm',
                    'e_24_mm',
                    'et_24_mm',
                    'et_ref_24_mm',
                    'se_root',
                    'biomass_prod',
                    'epoch_ends',
                    'epoch_starts']
        ds = PF.ds_remove_except(ds, keep_vars)
    elif isinstance(export_vars, list):
        keep_vars = copy.copy(export_vars)
        keep_vars = np.unique(keep_vars + ['epoch_ends', 'epoch_starts']).tolist()
        ds = PF.ds_remove_except(ds, keep_vars)
    else:
        raise ValueError

    ds = g.variables.fill_attrs(ds)
    ds = ds.transpose("epoch", "lat", "lon") # set dimension order the same for all vars.
    
    fp, fn = os.path.split(input_data)
    fn = fn.replace("_input", "_output")
    ds, fh = calculate_ds(ds, os.path.join(fp, fn), "--> Saving outputs.")

    if export_to_tif:
        files = PF.export_ds_to_tif(ds, keep_vars, None)
        ds.close()
        os.remove(fh)
        return files
    else:
        return ds

if __name__ == "__main__":

    # project_folder = r"/Volumes/Data/FAO/WaPOR_vs_pyWaPOR/pyWAPOR_v1"
    # startdate = date = "2021-07-01"

    # level = "level_1"
    et_look_version = "v2"
    output = None
    # input_data = None

    input_data = r"/Users/hmcoerver/pywapor_notebooks/level_1/et_look_input.nc"

    ds = main(input_data, et_look_version=et_look_version)
