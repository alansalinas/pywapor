#%%
import os
from osgeo import gdal
import numpy as np
import datetime
import warnings
import pywapor.et_look_dev as ETLook_dev
import pywapor.et_look_v2 as ETLook_v2
import pywapor.general.outputs as out
import pywapor.general.variables as vars
import pywapor.general.processing_functions as PF
import copy
import xarray as xr

def get_geoinfo(example_filepath):
    ds = gdal.Open(example_filepath)
    geo_ex = ds.GetGeoTransform()
    proj_ex = ds.GetProjection()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    dlat, dlon = PF.calc_dlat_dlon(geo_ex, xsize, ysize)
    dem_resolution = (np.nanmean(dlon) + np.nanmean(dlat))/2
    return dem_resolution, geo_ex, proj_ex

def open_array(key, value, date, folders):
    fp = create_fp(key, value, date, folders)
    if os.path.exists(fp):
        # array = create_array(fp)
        ds = gdal.Open(fp)
        band = ds.GetRasterBand(1)
        ndv = band.GetNoDataValue()
        array = band.ReadAsArray()
        array[array == ndv] = np.nan
        return array
    else:
        # print("'{0}' not found.".format(key))
        return None

def create_fp(key, value, date, folders):
    if value["time"] == "daily":
        date = date.strftime("_%Y%m%d")
    elif value["time"] == "static":
        date = ""
    elif value["time"] == "yearly":
        date = date.strftime("_%Y")
    else:
        print("ERROR: invalid value")

    fn = "{0}{1}.tif".format(key, date)
    fp = os.path.join(folders[value["time"]], fn)

    return fp

def main(project_folder, date, level = "level_1", et_look_version = "v2", 
            output = None, input_data = None):

    # Version
    if et_look_version == "v2":
        ETLook = ETLook_v2
        print("--> Running ETLook_v2")
    elif et_look_version == "dev":
        ETLook = ETLook_dev
        print("--> Running ETLOOK_dev")
    c = ETLook.constants

    # Date
    if isinstance(date, str):
        date = datetime.datetime.strptime(date, "%Y-%m-%d")
    date_str = date.strftime("%Y%m%d")

    # Disable Warnings
    warnings.filterwarnings('ignore')

    # Folders
    input_folder_date = os.path.join(project_folder, level, date_str)
    input_folder_static = os.path.join(project_folder, level, "static")
    folders = {"daily": input_folder_date,
               "static": input_folder_static,
               "yearly": input_folder_static,}
    output_folder_date = os.path.join(project_folder, f"out_{level}", date_str)
    if not os.path.exists(output_folder_date) and isinstance(input_data, type(None)):
        os.makedirs(output_folder_date)

    if isinstance(output, type(None)):
        out_final = out.outputs
    elif output == "all":
        out_final = {k: True for k in out.outputs.keys()}
    elif isinstance(output, dict):
        for key, value in output.items():
            out.outputs[key] = value
            out_final = out.outputs

    for key, value in vars.outputs.items():
        value["output"] = out_final[key]
        value["file_path"] = os.path.join(output_folder_date, 
        "{0}_{1}.tif".format(value["file_name"], date_str))

    # Inputs
    if isinstance(input_data, type(None)):
        id = {v["array_name"]: open_array(k, v, date, folders) for k, v in vars.inputs.items()}
    else:
        id = copy.deepcopy(input_data)

    spatial_or_constant = ["lw_offset", "lw_slope", "r0_bare", "r0_full", 
                            "rn_offset", "rn_slope", "t_opt", "vpd_slope"]

    # Logging
    for name, value in id.items():
        if name == "se_root" and isinstance(value, type(None)):
            print("----> se_root will be calculated from instantaneous data.")
        if name == "t_air_24" and isinstance(value, type(None)):
            print("----> t_air_24 will be calculated from t_air_min_24 and t_air_max_24.")
        if name == "ra_24" and isinstance(value, type(None)):
            print("----> ra_24 will be calculated from trans_24.")
        if name == "trans_24" and isinstance(value, type(None)):
            print("----> trans_24 will be calculated from ra_24.")
        if name == "z_oro" and isinstance(value, type(None)):
            print("----> z_oro will be calculated from slope.")
        if name == "u_24" and isinstance(value, type(None)):
            print("----> u_24 will be calculated from u2m_24 and v2m_24.")
        if name in spatial_or_constant and isinstance(value, type(None)):
            print(f"----> {name} will be constant at {getattr(c, name):.4f}.")

    # Resolution
    if isinstance(input_data, type(None)):
        example_filepath = create_fp("r0", vars.inputs["r0"], date, folders)
        geoinfo = get_geoinfo(example_filepath)
        resolution, geo_ex, proj_ex = geoinfo
        print(f"----> resolution is ~{resolution:.0f} meter.")

    for param in spatial_or_constant:
        if isinstance(id[param], type(None)):
            id[param] = np.ones_like(id["ndvi"]) * getattr(c, param)

    if isinstance(id["u_24"], type(None)):
        id["u_24"] = np.sqrt(id["v2m_24"]**2 + id["u2m_24"]**2)

    # Outputs
    od = dict()

    # Constants
    doy = int(date.strftime("%j"))
    sc = ETLook.solar_radiation.seasonal_correction(doy)
    decl = ETLook.solar_radiation.declination(doy)
    iesd = ETLook.solar_radiation.inverse_earth_sun_distance(doy)

    ######################## MODEL ETLOOK ####################################

    # **effective_leaf_area_index*********************************************
    od["vc"] = ETLook.leaf.vegetation_cover(id["ndvi"], c.nd_min, c.nd_max, c.vc_pow)
    od["lai"] = ETLook.leaf.leaf_area_index(od["vc"], c.vc_min, c.vc_max, c.lai_pow)
    od["lai_eff"] = ETLook.leaf.effective_leaf_area_index(od["lai"])
    
    # *******TRANSPIRATION COMPONENT******************************************
    
    # **soil fraction*********************************************************
    od["sf_soil"] = ETLook.radiation.soil_fraction(od["lai"])

    # **atmospheric canopy resistance***********************************************
    od["lat"] = ETLook.solar_radiation.latitude_rad(id["lat_deg"])
    od["ws"] = ETLook.solar_radiation.sunset_hour_angle(od["lat"], decl)
    if isinstance(id["ra_24"], type(None)):
        od["slope"] = ETLook.solar_radiation.slope_rad(id["slope_deg"])
        od["aspect"] = ETLook.solar_radiation.aspect_rad(id["aspect_deg"])
        od["ra_24_toa"] = ETLook.solar_radiation.daily_solar_radiation_toa_new(sc, decl, iesd, od["lat"], doy, od["slope"], od["aspect"])
        od["ra_24_toa_flat"] = ETLook.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, od["lat"], od["ws"])
        od["diffusion_index"] = ETLook.solar_radiation.diffusion_index(id["trans_24"], c.diffusion_slope, c.diffusion_intercept)
        od["ra_24"] = ETLook.solar_radiation.daily_total_solar_radiation(od["ra_24_toa"], od["ra_24_toa_flat"], od["diffusion_index"], id["trans_24"])
    else:
        od["ra_24"] = np.copy(id["ra_24"])
        # od["ra_24_toa"] = ETLook.solar_radiation.daily_solar_radiation_toa_new(sc, decl, iesd, od["lat"], doy, od["slope"], od["aspect"])
        od["ra_24_toa_flat"] = ETLook.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, od["lat"], od["ws"])
        od["trans_24"] = od["ra_24"] / od["ra_24_toa_flat"]
        id["trans_24"] = np.copy(od["trans_24"])

    od["stress_rad"] = ETLook.stress.stress_radiation(od["ra_24"])
    id["p_air_0_24"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_0_24"])
    od["p_air_24"] = ETLook.meteo.air_pressure_daily(id["z"], id["p_air_0_24"])
    od["vp_24"] = ETLook.meteo.vapour_pressure_from_specific_humidity_daily(id["qv_24"], od["p_air_24"])

    if isinstance(id["t_air_24"], type(None)):
        id["t_air_24"] = np.mean([id["t_air_min_24"], id["t_air_max_24"]], axis = 0)

    if et_look_version == "v2":
        od["svp_24"] = ETLook.meteo.saturated_vapour_pressure(id["t_air_24"])
    elif et_look_version == "dev":
        od["svp_24"] = ETLook.meteo.saturated_vapour_pressure_average(
                    ETLook.meteo.saturated_vapour_pressure_maximum(id["t_air_max_24"]),
                    ETLook.meteo.saturated_vapour_pressure_minimum(id["t_air_min_24"]))

    od["vpd_24"] = ETLook.meteo.vapour_pressure_deficit_daily(od["svp_24"], od["vp_24"])
    od["stress_vpd"] = ETLook.stress.stress_vpd(od["vpd_24"], id["vpd_slope"])
    od["stress_temp"] = ETLook.stress.stress_temperature(id["t_air_24"], id["t_opt"], c.t_min, c.t_max)

    if et_look_version == "dev":
        id["rs_min"] = np.where(id["land_mask"] == 3, 400, 100)

    od["r_canopy_0"] = ETLook.resistance.atmospheric_canopy_resistance(od["lai_eff"], od["stress_rad"], od["stress_vpd"], od["stress_temp"], id["rs_min"], c.rcan_max)

    # **net radiation canopy******************************************************
    od["t_air_k_24"] = ETLook.meteo.air_temperature_kelvin_daily(id["t_air_24"])
    od["l_net"] = ETLook.radiation.longwave_radiation_fao(od["t_air_k_24"], od["vp_24"], id["trans_24"], c.vp_slope, c.vp_offset, id["lw_slope"], id["lw_offset"])
    od["int_mm"] = ETLook.evapotranspiration.interception_mm(id["P_24"], od["vc"], od["lai"], c.int_max)
    od["lh_24"] = ETLook.meteo.latent_heat_daily(id["t_air_24"])
    od["int_wm2"] = ETLook.radiation.interception_wm2(od["int_mm"], od["lh_24"])
    od["rn_24"] = ETLook.radiation.net_radiation(id["r0"], od["ra_24"], od["l_net"], od["int_wm2"])
    od["rn_24_canopy"] = ETLook.radiation.net_radiation_canopy(od["rn_24"], od["sf_soil"])

    # **Soil Moisture Content***********************************************************
    if isinstance(id["se_root"], type(None)):
        if et_look_version == "dev":
            c.z0m_full = 0.04 + 0.01 * (resolution - 30)/(250-30)
            od["lst_zone_mean"] = lst_zone_mean(id, geo_ex, proj_ex, date, folders)
        od, id = se_root(id, od, ETLook, date)
    else:
        od["se_root"] = np.copy(id["se_root"])

    # find water region using ndvi
    if et_look_version == "dev":
        if resolution < 250:
            id["land_mask"] = np.where(id["land_mask"] == 2, 1, id["land_mask"])
            id["land_mask"] = np.where(id["ndvi"] < 0, 2, id["land_mask"])

    od["stress_moist"] = ETLook.stress.stress_moisture(od["se_root"], c.tenacity)
    od["r_canopy"] = ETLook.resistance.canopy_resistance(od["r_canopy_0"], od["stress_moist"], c.rcan_max)

    # **initial canopy aerodynamic resistance***********************************************************
    
    if isinstance(id["z_oro"], type(None)):
        od["slope"] = ETLook.solar_radiation.slope_rad(id["slope_deg"])
        od["z_oro"] = ETLook.roughness.orographic_roughness(od["slope"], resolution)
    else:
        od["z_oro"] = np.copy(id["z_oro"])

    od["z_obst"] = ETLook.roughness.obstacle_height(id["ndvi"], id["z_obst_max"], c.ndvi_obs_min, c.ndvi_obs_max, c.obs_fr)
    od["z0m"] = ETLook.roughness.roughness_length(od["lai"], od["z_oro"], od["z_obst"], id["z_obst_max"], id["land_mask"])
    od["ra_canopy_init"] = ETLook.neutral.initial_canopy_aerodynamic_resistance(id["u_24"], od["z0m"], c.z_obs)

    # **windspeed blending height daily***********************************************************
    od["u_b_24"] = ETLook.meteo.wind_speed_blending_height_daily(id["u_24"], c.z_obs, c.z_b)

    # **ETLook.neutral.initial_daily_transpiration***********************************************************
    od["ad_dry_24"] = ETLook.meteo.dry_air_density_daily(od["p_air_24"], od["vp_24"], od["t_air_k_24"])
    od["ad_moist_24"] = ETLook.meteo.moist_air_density_daily(od["vp_24"], od["t_air_k_24"])
    od["ad_24"] = ETLook.meteo.air_density_daily(od["ad_dry_24"], od["ad_moist_24"])
    od["psy_24"] = ETLook.meteo.psychrometric_constant_daily(od["p_air_24"], od["lh_24"])
    od["ssvp_24"] = ETLook.meteo.slope_saturated_vapour_pressure_daily(id["t_air_24"])
    od["t_24_init"] = ETLook.neutral.initial_daily_transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"], od["ra_canopy_init"])

    # **ETLook.unstable.initial_sensible_heat_flux_canopy_daily***********************************************************
    od["h_canopy_24_init"] = ETLook.unstable.initial_sensible_heat_flux_canopy_daily(od["rn_24_canopy"], od["t_24_init"])

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    od["disp"] = ETLook.roughness.displacement_height(od["lai"], od["z_obst"], id["land_mask"], c.c1)
    od["u_star_24_init"] = ETLook.unstable.initial_friction_velocity_daily(od["u_b_24"], od["z0m"], od["disp"], c.z_b)

    # **ETLook.unstable.transpiration***********************************************************
    od["t_24"] = ETLook.unstable.transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"], od["h_canopy_24_init"], od["t_air_k_24"], od["u_star_24_init"], od["z0m"], od["disp"], od["u_b_24"], c.z_obs, c.z_b, c.iter_h)
    od["t_24_mm"] = ETLook.unstable.transpiration_mm(od["t_24"], od["lh_24"])

    if et_look_version == "dev":
        od["tpot_24"] = ETLook.unstable.transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"] * od["stress_moist"], od["h_canopy_24_init"], od["t_air_k_24"], od["u_star_24_init"], od["z0m"], od["disp"], od["u_b_24"], c.z_obs, c.z_b, c.iter_h)
        od["tpot_24_mm"] = ETLook.unstable.transpiration_mm(od["tpot_24"], od["lh_24"])

    #*******EVAPORATION COMPONENT****************************************************************

    # **ETLook.radiation.net_radiation_soil***********************************************************
    od["rn_24_soil"] = ETLook.radiation.net_radiation_soil(od["rn_24"], od["sf_soil"])

    # **ETLook.resistance.soil_resistance***********************************************************
    od["r_soil"] = ETLook.resistance.soil_resistance(od["se_root"], id["land_mask"], c.r_soil_pow, c.r_soil_min)

    # **ETLook.resistance.soil_resistance***********************************************************
    od["ra_soil_init"] = ETLook.neutral.initial_soil_aerodynamic_resistance(id["u_24"], c.z_obs)

    # **ETLook.unstable.initial_friction_velocity_soil_daily***********************************************************
    od["u_star_24_soil_init"] = ETLook.unstable.initial_friction_velocity_soil_daily(od["u_b_24"], od["disp"], c.z_b)

    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************
    if et_look_version == "dev":
        stc = ETLook.radiation.soil_thermal_conductivity(c.se_top)
        vhc = ETLook.radiation.volumetric_heat_capacity(c.se_top, c.porosity)
    elif et_look_version == "v2":
        stc = ETLook.radiation.soil_thermal_conductivity(od["se_root"])
        vhc = ETLook.radiation.volumetric_heat_capacity(od["se_root"], c.porosity)
       
    dd = ETLook.radiation.damping_depth(stc, vhc)
    od["g0_bs"] = ETLook.radiation.bare_soil_heat_flux(doy, dd, stc, id["t_amp_year"], od["lat"])
    od["g0_24"] = ETLook.radiation.soil_heat_flux(od["g0_bs"], od["sf_soil"], id["land_mask"], od["rn_24_soil"], id["trans_24"], od["ra_24"], od["l_net"], id["rn_slope"], id["rn_offset"])
    od["e_24_init"] = ETLook.neutral.initial_daily_evaporation(od["rn_24_soil"], od["g0_24"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_soil"], od["ra_soil_init"])
    od["h_soil_24_init"] = ETLook.unstable.initial_sensible_heat_flux_soil_daily(od["rn_24_soil"], od["e_24_init"], od["g0_24"])

    # **ETLook.unstable.evaporation***********************************************************
    od["e_24"] = ETLook.unstable.evaporation(od["rn_24_soil"], od["g0_24"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_soil"], od["h_soil_24_init"], od["t_air_k_24"], od["u_star_24_soil_init"], od["disp"], od["u_b_24"], c.z_b, c.z_obs, c.iter_h)
    od["e_24_mm"] = ETLook.unstable.evaporation_mm(od["e_24"], od["lh_24"])
    od["et_24_mm"] = ETLook.evapotranspiration.et_actual_mm(od["e_24_mm"], od["t_24_mm"])
    
    if et_look_version == "dev":
        od["e_24_mm"][od["e_24_mm"] < 0] = 0
        od["et_24_mm"][od["et_24_mm"] < 0] = 0

    # **ETLook.unstable.evaporation***********************************************************
    od["rn_24_grass"] = ETLook.radiation.net_radiation_grass(od["ra_24"], od["l_net"], c.r0_grass)
    od["et_ref_24"] = ETLook.evapotranspiration.et_reference(od["rn_24_grass"], od["ad_24"], od["psy_24"], od["vpd_24"], od["ssvp_24"], id["u_24"])
    od["et_ref_24_mm"] = ETLook.evapotranspiration.et_reference_mm(od["et_ref_24"], od["lh_24"])

    if et_look_version == "v2":
        od["et_ref_24_mm"][od["et_ref_24_mm"] < 0] = 0

    if et_look_version == "dev":
        eps_a = ETLook.stress.epsilon_autotrophic_respiration()     
        od["lue"] = ETLook.biomass.lue(id["lue_max"], od["stress_temp"], od["stress_moist"], eps_a)
        od["fpar"] = ETLook.leaf.fpar(od["vc"], id["ndvi"])
        od["apar"] = ETLook.leaf.apar(od["ra_24"], od["fpar"])       
        od["biomass_prod"] = ETLook.biomass.biomass(od["apar"], od["lue"])     

    if isinstance(input_data, type(None)):
        all_files = dict()
        example_filepath = create_fp("r0", vars.inputs["r0"], date, folders)
        geo_ex, proj_ex = get_geoinfo(example_filepath)[1:3]
        for var in vars.outputs.keys():
            if var in od.keys():
                if vars.outputs[var]["output"]:
                    PF.Save_as_tiff(vars.outputs[var]["file_path"], od[var], geo_ex, proj_ex)
                    all_files[var] = vars.outputs[var]["file_path"]
        return all_files
    else:
        return od, id

def se_root(id, od, ETLook, date, verbose = False):

    et_look_version = ETLook.__name__.split(".")[-1].split("_")[-1]

    # Version
    if et_look_version == "v2":
        # ETLook = ETLook_v2
        if not verbose:
            print("--> Running SEroot_v2")
    elif et_look_version == "dev":
        # ETLook = ETLook_dev
        if not verbose:
            print("--> Running SEroot_dev")
    c = ETLook.constants

    if not isinstance(id, dict):
        od = id

    # resolution, geo_ex, proj_ex = geoinfo

    doy = int(date.strftime("%j"))
    sc = ETLook.solar_radiation.seasonal_correction(doy)
    decl = ETLook.solar_radiation.declination(doy)
    day_angle = ETLook.clear_sky_radiation.day_angle(doy)
    dtime = date.hour + (date.minute / 60)

    if "r0_bare" not in id.keys():
        if isinstance(id, xr.Dataset):
            id["r0_bare"] = xr.ones_like(id["ndvi"]) * c.r0_bare
        else:
            id["r0_bare"] = np.ones_like(id["ndvi"]) * c.r0_bare
    if "r0_full" not in id.keys():
        if isinstance(id, xr.Dataset):
            id["r0_full"] = xr.ones_like(id["ndvi"]) * c.r0_full
        else:
            id["r0_full"] = np.ones_like(id["ndvi"]) * c.r0_full

    if et_look_version == "dev":
        id["t_air_i"][id["t_air_i"] < -270] = np.nan
    id["p_air_i"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_i"])
    id["p_air_0_i"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_0_i"])

    # if et_look_version == "dev":
    #     z0m_full = 0.04 + 0.01 * (resolution - 30)/(250-30)
    # elif et_look_version == "v2":
    #     z0m_full = 0.1

    od["vc"] = ETLook.leaf.vegetation_cover(id["ndvi"], c.nd_min, c.nd_max, c.vc_pow)

    od["t_air_k_i"] = ETLook.meteo.air_temperature_kelvin_inst(id["t_air_i"])
    od["vp_i"] = ETLook.meteo.vapour_pressure_from_specific_humidity_inst(id["qv_i"], id["p_air_i"])
    od["ad_moist_i"] = ETLook.meteo.moist_air_density_inst(od["vp_i"], od["t_air_k_i"])
    od["ad_dry_i"] = ETLook.meteo.dry_air_density_inst(id["p_air_i"], od["vp_i"], od["t_air_k_i"])
    od["ad_i"] = ETLook.meteo.air_density_inst(od["ad_dry_i"], od["ad_moist_i"])
    od["u_b_i_bare"] = ETLook.soil_moisture.wind_speed_blending_height_bare(id["u_i"], c.z0m_bare, c.z_obs, c.z_b)
    od["lon"] = ETLook.solar_radiation.longitude_rad(id["lon_deg"])
    od["lat"] = ETLook.solar_radiation.latitude_rad(id["lat_deg"])
    od["ha"] = ETLook.solar_radiation.hour_angle(sc, dtime, od["lon"])
    
    I0 = ETLook.clear_sky_radiation.solar_constant()
    ied = ETLook.clear_sky_radiation.inverse_earth_sun_distance(day_angle)
    od["h0"] = ETLook.clear_sky_radiation.solar_elevation_angle(od["lat"], decl, od["ha"])
    od["h0ref"] = ETLook.clear_sky_radiation.solar_elevation_angle_refracted(od["h0"])
    od["m"] = ETLook.clear_sky_radiation.relative_optical_airmass(id["p_air_i"], id["p_air_0_i"], od["h0ref"])
    od["rotm"] = ETLook.clear_sky_radiation.rayleigh_optical_thickness(od["m"])
    od["Tl2"] = ETLook.clear_sky_radiation.linke_turbidity(id["wv_i"], c.aod550_i, id["p_air_i"], id["p_air_0_i"])
    G0 = ETLook.clear_sky_radiation.extraterrestrial_irradiance_normal(I0, ied)
    od["B0c"] = ETLook.clear_sky_radiation.beam_irradiance_normal_clear(G0, od["Tl2"], od["m"], od["rotm"], od["h0"])
    od["Bhc"] = ETLook.clear_sky_radiation.beam_irradiance_horizontal_clear(od["B0c"], od["h0"])
    od["Dhc"] = ETLook.clear_sky_radiation.diffuse_irradiance_horizontal_clear(G0, od["Tl2"], od["h0"])
    
    od["ra_hor_clear_i"] = ETLook.clear_sky_radiation.ra_clear_horizontal(od["Bhc"], od["Dhc"])
    od["emiss_atm_i"] = ETLook.soil_moisture.atmospheric_emissivity_inst(od["vp_i"], od["t_air_k_i"])
    
    od["rn_bare"] = ETLook.soil_moisture.net_radiation_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], id["lst"], id["r0_bare"])
    od["rn_full"] = ETLook.soil_moisture.net_radiation_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], id["lst"], id["r0_full"])
    od["h_bare"] = ETLook.soil_moisture.sensible_heat_flux_bare(od["rn_bare"], c.fraction_h_bare)
    od["h_full"] = ETLook.soil_moisture.sensible_heat_flux_full(od["rn_full"], c.fraction_h_full)
    od["u_b_i_full"] = ETLook.soil_moisture.wind_speed_blending_height_full_inst(id["u_i"], c.z0m_full, c.z_obs, c.z_b)
    
    od["u_star_i_bare"] = ETLook.soil_moisture.friction_velocity_bare_inst(od["u_b_i_bare"], c.z0m_bare, c.disp_bare, c.z_b)
    od["u_star_i_full"] = ETLook.soil_moisture.friction_velocity_full_inst(od["u_b_i_full"], c.z0m_full, c.disp_full, c.z_b)
    od["L_bare"] = ETLook.soil_moisture.monin_obukhov_length_bare(od["h_bare"], od["ad_i"], od["u_star_i_bare"], od["t_air_k_i"])
    od["L_full"] = ETLook.soil_moisture.monin_obukhov_length_full(od["h_full"], od["ad_i"], od["u_star_i_full"], od["t_air_k_i"])
    
    od["u_i_soil"] = ETLook.soil_moisture.wind_speed_soil_inst(id["u_i"], od["L_bare"], c.z_obs)
    od["ras"] = ETLook.soil_moisture.aerodynamical_resistance_soil(od["u_i_soil"])
    od["raa"] = ETLook.soil_moisture.aerodynamical_resistance_bare(id["u_i"], od["L_bare"], c.z0m_bare, c.disp_bare, c.z_obs)
    od["rac"] = ETLook.soil_moisture.aerodynamical_resistance_full(id["u_i"], od["L_full"], c.z0m_full, c.disp_full, c.z_obs)
    
    od["t_max_bare"] = ETLook.soil_moisture.maximum_temperature_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["raa"], od["ras"], id["r0_bare"])
    od["t_max_full"] = ETLook.soil_moisture.maximum_temperature_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["rac"], id["r0_full"])

    # od["w_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])
    # od["t_dew_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])
    # od["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst(id["t_air_i"], od["t_dew_i"], id["p_air_i"])
    od["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst_new(id["t_air_i"], id["qv_i"], id["p_air_i"])
    #t_wet_i_new2 = ETLook.soil_moisture.wet_bulb_temperature_inst_new2(t_air_i)

    od["lst_max"] = ETLook.soil_moisture.maximum_temperature(od["t_max_bare"], od["t_max_full"], od["vc"])

    if et_look_version == "v2":
        od["t_wet_k_i"] = ETLook.meteo.wet_bulb_temperature_kelvin_inst(od["t_wet_i"])
        od["lst_min"] = ETLook.soil_moisture.minimum_temperature(od["t_wet_k_i"], od["t_air_k_i"], od["vc"])
    elif et_look_version == "dev":
        od["t_min_bare"] = ETLook.soil_moisture.minimum_temperature_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["raa"], od["ras"], od["lst_zone_mean"], c.r0_bare_wet)
        od["t_min_full"] = ETLook.soil_moisture.minimum_temperature_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["rac"], od["lst_zone_mean"], id["r0_full"])  
        od["lst_min"] = ETLook.soil_moisture.maximum_temperature(od["t_min_bare"], od["t_min_full"], od["vc"])

    od["se_root"] = ETLook.soil_moisture.soil_moisture_from_maximum_temperature(od["lst_max"], id["lst"], od["lst_min"])

    return od, id

def lst_zone_mean(id, geo_ex, proj_ex, date, folders):
    size_y, size_x = id["ndvi"].shape
    size_y_zone = int(np.ceil(size_y/200))
    size_x_zone = int(np.ceil(size_x/200)) 
    array_fake = np.ones([size_y_zone, size_x_zone])
    geo_new = tuple([geo_ex[0], geo_ex[1] * 200, geo_ex[2], geo_ex[3], geo_ex[4], geo_ex[5]*200])
    MEM_file = PF.Save_as_MEM(array_fake, geo_new, proj_ex)

    LST_filename = create_fp("LST", vars.inputs["LST"], date, folders)
    
    dest_lst_zone_large = PF.reproject_dataset_example(LST_filename, MEM_file, 4)
    lst_zone_mean_large = dest_lst_zone_large.GetRasterBand(1).ReadAsArray()
    lst_zone_mean_large[lst_zone_mean_large==0] = -9999
    lst_zone_mean_large[np.isnan(lst_zone_mean_large)] = -9999

    if np.nanmax(lst_zone_mean_large) == -9999:
        for x in range(0, size_x_zone):
            for y in range(0, size_y_zone):
                lst_zone_mean_large[y, x] = np.nanmean(id["lst"][y*200:np.minimum((y+1)*200, size_y-1), x*200:np.minimum((x+1)*200, size_x-1)])
        lst_zone_mean_large[np.isnan(lst_zone_mean_large)] = -9999

    lst_zone_mean_large = PF.gap_filling(lst_zone_mean_large, -9999, 1)
    dest_lst_zone_large = PF.Save_as_MEM(lst_zone_mean_large, geo_new, proj_ex)
    
    dest_lst_zone = PF.reproject_dataset_example(dest_lst_zone_large, LST_filename, 6)
    lst_zone_mean = dest_lst_zone.GetRasterBand(1).ReadAsArray()
    return lst_zone_mean

# if __name__ == "__main__":

#     project_folder = r"/Volumes/Data/FAO/WaPOR_vs_pyWaPOR/pyWAPOR_v1"
#     startdate = date = "2021-07-01"

#     level = "level_1"
#     et_look_version = "v2"
#     output = None
#     input_data = None

#     # all_files = main(project_folder, startdate)
