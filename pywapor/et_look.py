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

def main(project_folder, date, level = "level_1", et_look_version = "v2"):

    if et_look_version == "v2":
        ETLook = ETLook_v2
    elif et_look_version == "dev":
        ETLook = ETLook_dev
    c = ETLook.constants

    if isinstance(date, str):
        date = datetime.datetime.strptime(date, "%Y-%m-%d")

    # Do not show warnings
    warnings.filterwarnings('ignore')

    # Define date string
    date_str = "%d%02d%02d" %(date.year, date.month, date.day)

    # Input folder date
    input_folder_date = os.path.join(project_folder, level, date_str)
    input_folder_static = os.path.join(project_folder, level, "Static")

    ############################ Define inputs ################################

    folders = {"daily": input_folder_date,
               "static": input_folder_static,
               "yearly": input_folder_static,}

    def create_array(fp):
        ds = gdal.Open(fp)
        return ds.GetRasterBand(1).ReadAsArray()

    def create_fp(key, value, date, folders = folders):
        if value["time"] == "daily":
            date = "_" + date_str
        elif value["time"] == "static":
            date = ""
        elif value["time"] == "yearly":
            date = "_" + str(date.year)
        else:
            print("ERROR: invalid value")

        fn = "{0}{1}.tif".format(key, date)
        fp = os.path.join(folders[value["time"]], fn)

        if not os.path.isfile(fp):
            print(f"WARNING: {fp} does not exists")

        return fp

    # Create QC array
    fp = create_fp("LST", vars.inputs["LST"], date)
    # fp = create_fp("NDVI", vars.inputs["NDVI"])
    init = create_array(fp)
    QC = np.ones(init.shape)
    QC[init == -9999] = np.nan

    def open_array(key, value, date, folders = folders, mask = QC):
        fp = create_fp(key, value, date, folders = folders)
        if os.path.exists(fp):
            array = create_array(fp)
            array[np.isnan(mask)] = np.nan
            return array
        else:
            print("'{0}' not found.".format(key))
            return None

    # Input Data
    id = {v["array_name"]: open_array(k, v, date) for k, v in vars.inputs.items()}

    ############################ Define outputs ###############################

    # Output folder date
    output_folder_date = os.path.join(project_folder, f"out_{level}", date_str)
    if not os.path.exists(output_folder_date):
        os.makedirs(output_folder_date)

    if isinstance(id["se_root"], type(None)):
        method = 1
    else:
        method = 2

    for key, value in vars.outputs.items():
        value["output"] = bool(getattr(out, key))
        value["file_path"] = os.path.join(output_folder_date, 
        "{0}_{1}.tif".format(value["file_name"], date_str))

    od = dict()

    if "se_root" in id.keys():
        od["se_root"] = id["se_root"]

    id["p_air_0_24"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_0_24"])

    # example file
    
    def get_geoinfo(example_filepath):
        ds = gdal.Open(example_filepath)
        geo_ex = ds.GetGeoTransform()
        proj_ex = ds.GetProjection()
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        dlat, dlon = PF.calc_dlat_dlon(geo_ex, xsize, ysize)
        dem_resolution = (np.nanmean(dlon) + np.nanmean(dlat))/2
        return dem_resolution, geo_ex, proj_ex

    if method == 1:

        if et_look_version == "dev":
            id["t_air_i"][id["t_air_i"] < -270] = np.nan
        id["p_air_i"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_i"])
        id["p_air_0_i"] = ETLook.meteo.air_pressure_kpa2mbar(id["p_air_0_i"])

    doy = int(date.strftime("%j"))

    if et_look_version == "dev":
        example_filepath = create_fp("ALBEDO", vars.inputs["ALBEDO"], date)
        dem_resolution, geo_ex, proj_ex = get_geoinfo(example_filepath)
        id["rs_min"] = np.where(id["land_mask"] == 3, 400, 100)
        z0m_full = 0.04 + 0.01 * (dem_resolution - 30)/(250-30)
    elif et_look_version == "v2":
        z0m_full = 0.1

    ######################## MODEL ETLOOK #########################################

    # **effective_leaf_area_index**************************************************
    od["vc"] = ETLook.leaf.vegetation_cover(id["ndvi"], c.nd_min, c.nd_max, c.vc_pow)
    od["lai"] = ETLook.leaf.leaf_area_index(od["vc"], c.vc_min, c.vc_max, c.lai_pow)
    od["lai_eff"] = ETLook.leaf.effective_leaf_area_index(od["lai"])
    #*******TRANSPIRATION COMPONENT****************************************************************
    # **soil fraction************************************************************** 
    od["sf_soil"] = ETLook.radiation.soil_fraction(od["lai"])

    # **atmospheric canopy resistance***********************************************
    iesd = ETLook.solar_radiation.inverse_earth_sun_distance(doy)

    sc = ETLook.solar_radiation.seasonal_correction(doy)

    b1 = 2 * np.pi * (doy - 81) / 365.0
    b2 = 2 * np.pi * (doy - 81) / 364.0
    sc_dev = 0.1645 * np.sin(2 * b1) - 0.1255 * np.cos(b1) - 0.025 * np.sin(b1)
    sc_v2 =  0.1645 * np.sin(2 * b2) - 0.1255 * np.cos(b2) - 0.025 * np.sin(b2)
    if sc == sc_dev:
        print("USING ETLOOK_dev")
        assert et_look_version == "dev", "using wrong ETLook import."
    elif sc == sc_v2:
        print("USING ETLook_v2")
        assert et_look_version == "v2", "using wrong ETLook import."

    day_angle = ETLook.clear_sky_radiation.day_angle(doy)
    
    decl = ETLook.solar_radiation.declination(doy)

    od["lat"] = ETLook.solar_radiation.latitude_rad(id["lat_deg"])
    od["slope"] = ETLook.solar_radiation.slope_rad(id["slope_deg"])
    od["aspect"] = ETLook.solar_radiation.aspect_rad(id["aspect_deg"])

    #ra_24_toa = ETLook.solar_radiation.daily_solar_radiation_toa(sc, decl, iesd, lat, slope, aspect)
    od["ra_24_toa"] = ETLook.solar_radiation.daily_solar_radiation_toa_new(sc, decl, iesd, od["lat"], doy, od["slope"], od["aspect"])
    
    od["ws"] = ETLook.solar_radiation.sunset_hour_angle(od["lat"], decl)
    od["ra_24_toa_flat"] = ETLook.solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, od["lat"], od["ws"])
    od["diffusion_index"] = ETLook.solar_radiation.diffusion_index(id["trans_24"], c.diffusion_slope, c.diffusion_intercept)

    # choose one of the two options below
    #ra_24 = ETLook.solar_radiation.daily_solar_radiation_flat(ra_24_toa_flat, trans_24)
    od["ra_24"] = ETLook.solar_radiation.daily_total_solar_radiation(od["ra_24_toa"], od["ra_24_toa_flat"], od["diffusion_index"], id["trans_24"])
    od["stress_rad"] = ETLook.stress.stress_radiation(od["ra_24"])
    od["p_air_24"] = ETLook.meteo.air_pressure_daily(id["z"], id["p_air_0_24"])
    od["vp_24"] = ETLook.meteo.vapour_pressure_from_specific_humidity_daily(id["qv_24"], od["p_air_24"])
    od["svp_24"] = ETLook.meteo.saturated_vapour_pressure_average(
                ETLook.meteo.saturated_vapour_pressure_maximum(id["t_air_max_24"]),
                ETLook.meteo.saturated_vapour_pressure_minimum(id["t_air_min_24"]))
    od["vpd_24"] = ETLook.meteo.vapour_pressure_deficit_daily(od["svp_24"], od["vp_24"])
    od["stress_vpd"] = ETLook.stress.stress_vpd(od["vpd_24"], c.vpd_slope)
    od["stress_temp"] = ETLook.stress.stress_temperature(id["t_air_24"], c.t_opt, c.t_min, c.t_max)
    od["r_canopy_0"] = ETLook.resistance.atmospheric_canopy_resistance(od["lai_eff"], od["stress_rad"], od["stress_vpd"], od["stress_temp"], id["rs_min"], c.rcan_max)

    # del od["aspect"], od["ra_24_toa"], od["ws"], od["ra_24_toa_flat"], od["diffusion_index"], od["svp_24"]

    # **net radiation canopy******************************************************
    od["t_air_k_24"] = ETLook.meteo.air_temperature_kelvin_daily(id["t_air_24"])
    # select one of the below two
    #l_net = ETLook.radiation.longwave_radiation_fao_etref(t_air_k_24, vp_24, trans_24)
    od["l_net"] = ETLook.radiation.longwave_radiation_fao(od["t_air_k_24"], od["vp_24"], id["trans_24"], c.vp_slope, c.vp_offset, c.lw_slope, c.lw_offset)
    od["int_mm"] = ETLook.evapotranspiration.interception_mm(id["P_24"], od["vc"], od["lai"], c.int_max)
    od["lh_24"] = ETLook.meteo.latent_heat_daily(id["t_air_24"])
    od["int_wm2"] = ETLook.radiation.interception_wm2(od["int_mm"], od["lh_24"])
    od["rn_24"] = ETLook.radiation.net_radiation(id["r0"], od["ra_24"], od["l_net"], od["int_wm2"])
    od["rn_24_canopy"] = ETLook.radiation.net_radiation_canopy(od["rn_24"], od["sf_soil"])

    # del od["int_mm"], od["int_wm2"]

    # **canopy resistance***********************************************************

    if method == 1:
        
        od["t_air_k_i"] = ETLook.meteo.air_temperature_kelvin_inst(id["t_air_i"])
        od["vp_i"] = ETLook.meteo.vapour_pressure_from_specific_humidity_inst(id["qv_i"], id["p_air_i"])
        od["ad_moist_i"] = ETLook.meteo.moist_air_density_inst(od["vp_i"], od["t_air_k_i"])
        od["ad_dry_i"] = ETLook.meteo.dry_air_density_inst(id["p_air_i"], od["vp_i"], od["t_air_k_i"])
        od["ad_i"] = ETLook.meteo.air_density_inst(od["ad_dry_i"], od["ad_moist_i"])
        od["u_b_i_bare"] = ETLook.soil_moisture.wind_speed_blending_height_bare(id["u_i"], c.z0m_bare, c.z_obs, c.z_b)
        od["lon"] = ETLook.solar_radiation.longitude_rad(id["lon_deg"])
        od["ha"] = ETLook.solar_radiation.hour_angle(sc, id["dtime"], od["lon"])
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
        
        od["rn_bare"] = ETLook.soil_moisture.net_radiation_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], id["lst"], c.r0_bare)
        od["rn_full"] = ETLook.soil_moisture.net_radiation_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], id["lst"], c.r0_full)
        od["h_bare"] = ETLook.soil_moisture.sensible_heat_flux_bare(od["rn_bare"], c.fraction_h_bare)
        od["h_full"] = ETLook.soil_moisture.sensible_heat_flux_full(od["rn_full"], c.fraction_h_full)
        od["u_b_i_full"] = ETLook.soil_moisture.wind_speed_blending_height_full_inst(id["u_i"], z0m_full, c.z_obs, c.z_b)
        
        od["u_star_i_bare"] = ETLook.soil_moisture.friction_velocity_bare_inst(od["u_b_i_bare"], c.z0m_bare, c.disp_bare, c.z_b)
        od["u_star_i_full"] = ETLook.soil_moisture.friction_velocity_full_inst(od["u_b_i_full"], z0m_full, c.disp_full, c.z_b)
        od["L_bare"] = ETLook.soil_moisture.monin_obukhov_length_bare(od["h_bare"], od["ad_i"], od["u_star_i_bare"], od["t_air_k_i"])
        od["L_full"] = ETLook.soil_moisture.monin_obukhov_length_full(od["h_full"], od["ad_i"], od["u_star_i_full"], od["t_air_k_i"])
        
        od["u_i_soil"] = ETLook.soil_moisture.wind_speed_soil_inst(id["u_i"], od["L_bare"], c.z_obs)
        od["ras"] = ETLook.soil_moisture.aerodynamical_resistance_soil(od["u_i_soil"])
        od["raa"] = ETLook.soil_moisture.aerodynamical_resistance_bare(id["u_i"], od["L_bare"], c.z0m_bare, c.disp_bare, c.z_obs)
        od["rac"] = ETLook.soil_moisture.aerodynamical_resistance_full(id["u_i"], od["L_full"], z0m_full, c.disp_full, c.z_obs)
        
        # Create MEM file zones
        if et_look_version == "dev":

            size_y, size_x = QC.shape
            size_y_zone = int(np.ceil(size_y/200))
            size_x_zone = int(np.ceil(size_x/200)) 
            array_fake = np.ones([size_y_zone, size_x_zone])
            geo_new = tuple([geo_ex[0], geo_ex[1] * 200, geo_ex[2], geo_ex[3], geo_ex[4], geo_ex[5]*200])
            MEM_file = PF.Save_as_MEM(array_fake, geo_new, proj_ex)

            LST_filename = create_fp("LST", vars.inputs["LST"], date)
            
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
            od["lst_zone_mean"] = dest_lst_zone.GetRasterBand(1).ReadAsArray()

        od["t_max_bare"] = ETLook.soil_moisture.maximum_temperature_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["raa"], od["ras"], c.r0_bare)
        od["t_max_full"] = ETLook.soil_moisture.maximum_temperature_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["rac"], c.r0_full)

        od["w_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])
        od["t_dew_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])

        #t_wet_i =      ETLook.soil_moisture.wet_bulb_temperature_inst        (t_air_i, t_dew_i, p_air_i)
        od["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst_new(id["t_air_i"], id["qv_i"], id["p_air_i"])
        #t_wet_i_new2 = ETLook.soil_moisture.wet_bulb_temperature_inst_new2(t_air_i)
   
        od["lst_max"] = ETLook.soil_moisture.maximum_temperature(od["t_max_bare"], od["t_max_full"], od["vc"])

        if et_look_version == "v2":
            od["t_wet_k_i"] = ETLook.meteo.wet_bulb_temperature_kelvin_inst(od["t_wet_i"])
            od["lst_min"] = ETLook.soil_moisture.minimum_temperature(od["t_wet_k_i"], od["t_air_k_i"], od["vc"])
        elif et_look_version == "dev":
            od["t_min_bare"] = ETLook.soil_moisture.minimum_temperature_bare(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["raa"], od["ras"], od["lst_zone_mean"], c.r0_bare_wet)
            od["t_min_full"] = ETLook.soil_moisture.minimum_temperature_full(od["ra_hor_clear_i"], od["emiss_atm_i"], od["t_air_k_i"], od["ad_i"], od["rac"], od["lst_zone_mean"], c.r0_full)  
            od["lst_min"] = ETLook.soil_moisture.maximum_temperature(od["t_min_bare"], od["t_min_full"], od["vc"])

        od["se_root"] = ETLook.soil_moisture.soil_moisture_from_maximum_temperature(od["lst_max"], id["lst"], od["lst_min"])

        # del od["t_air_k_i"], od["vp_i"], od["ad_moist_i"], od["ad_dry_i"], od["ad_i"], od["u_b_i_bare"]
        # del od["lon"], od["ha"], od["h0"], od["h0ref"], od["m"], od["rotm"], od["Tl2"], od["B0c"], od["Bhc"], od["Dhc"], od["ra_hor_clear_i"], od["emiss_atm_i"], od["rn_bare"], od["rn_full"], od["u_b_i_full"], od["u_star_i_bare"]
        # del od["u_star_i_full"], od["u_i_soil"], od["ras"], od["raa"], od["rac"], od["t_max_bare"], od["t_max_full"], od["w_i"], od["t_dew_i"], od["t_wet_i"], od["t_wet_k_i"], od["lst_max"]

    od["stress_moist"] = ETLook.stress.stress_moisture(od["se_root"], c.tenacity)
    od["r_canopy"] = ETLook.resistance.canopy_resistance(od["r_canopy_0"], od["stress_moist"], c.rcan_max)

    # del od["lai_eff"], od["stress_rad"], od["stress_vpd"]
    
    # **initial canopy aerodynamic resistance***********************************************************
    
    # find water region using ndvi
    if et_look_version == "dev":
        if dem_resolution < 250:
            id["land_mask"] = np.where(id["land_mask"] == 2, 1, id["land_mask"])
            id["land_mask"] = np.where(id["ndvi"] < 0, 2, id["land_mask"])
    elif et_look_version == "v2":
        example_filepath = create_fp("ALBEDO", vars.inputs["ALBEDO"], date)
        geo_ex, proj_ex = get_geoinfo(example_filepath)[1:3]
        dlat, dlon = PF.calc_dlat_dlon(geo_ex, id["ndvi"].shape[1], id["ndvi"].shape[0])
        dem_resolution = (np.nanmean(dlon) +np.nanmean(dlat))/2
    
    od["z_obst"] = ETLook.roughness.obstacle_height(id["ndvi"], id["z_obst_max"], c.ndvi_obs_min, c.ndvi_obs_max, c.obs_fr)
    od["z_oro"] = ETLook.roughness.orographic_roughness(od["slope"], dem_resolution)
    od["z0m"] = ETLook.roughness.roughness_length(od["lai"], od["z_oro"], od["z_obst"], id["z_obst_max"], id["land_mask"])
    od["ra_canopy_init"] = ETLook.neutral.initial_canopy_aerodynamic_resistance(id["u_24"], od["z0m"], c.z_obs)

    # del od["slope"], od["z_oro"]

    # **windspeed blending height daily***********************************************************
    od["u_b_24"] = ETLook.meteo.wind_speed_blending_height_daily(id["u_24"], c.z_obs, c.z_b)

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    od["disp"] = ETLook.roughness.displacement_height(od["lai"], od["z_obst"], id["land_mask"], c.c1)
    od["u_star_24_init"] = ETLook.unstable.initial_friction_velocity_daily(od["u_b_24"], od["z0m"], od["disp"], c.z_b)

    # del od["z_obst"]

    # **ETLook.neutral.initial_daily_transpiration***********************************************************
    od["ad_dry_24"] = ETLook.meteo.dry_air_density_daily(od["p_air_24"], od["vp_24"], od["t_air_k_24"])
    od["ad_moist_24"] = ETLook.meteo.moist_air_density_daily(od["vp_24"], od["t_air_k_24"])
    od["ad_24"] = ETLook.meteo.air_density_daily(od["ad_dry_24"], od["ad_moist_24"])
    od["psy_24"] = ETLook.meteo.psychrometric_constant_daily(od["p_air_24"], od["lh_24"])
    od["ssvp_24"] = ETLook.meteo.slope_saturated_vapour_pressure_daily(id["t_air_24"])
    od["t_24_init"] = ETLook.neutral.initial_daily_transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"], od["ra_canopy_init"])

    # del od["p_air_24"], od["vp_24"], od["ra_canopy_init"], od["ad_dry_24"], od["ad_moist_24"]

    # **ETLook.unstable.initial_sensible_heat_flux_canopy_daily***********************************************************
    od["h_canopy_24_init"] = ETLook.unstable.initial_sensible_heat_flux_canopy_daily(od["rn_24_canopy"], od["t_24_init"])

    # del od["t_24_init"]

    # **ETLook.unstable.transpiration***********************************************************

    od["t_24"] = ETLook.unstable.transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"], od["h_canopy_24_init"], od["t_air_k_24"], od["u_star_24_init"], od["z0m"], od["disp"], od["u_b_24"], c.z_obs, c.z_b, c.iter_h)
    od["t_24_mm"] = ETLook.unstable.transpiration_mm(od["t_24"], od["lh_24"])

    if et_look_version == "dev":
        od["tpot_24"] = ETLook.unstable.transpiration(od["rn_24_canopy"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_canopy"] * od["stress_moist"], od["h_canopy_24_init"], od["t_air_k_24"], od["u_star_24_init"], od["z0m"], od["disp"], od["u_b_24"], c.z_obs, c.z_b, c.iter_h)
        od["tpot_24_mm"] = ETLook.unstable.transpiration_mm(od["tpot_24"], od["lh_24"])

    # del od["rn_24_canopy"], od["r_canopy"], od["z0m"], od["u_star_24_init"], od["h_canopy_24_init"], od["t_24"]

    #*******EVAPORATION COMPONENT****************************************************************

    # **ETLook.radiation.net_radiation_soil***********************************************************
    od["sf_soil"] = ETLook.radiation.soil_fraction(od["lai"])
    od["rn_24_soil"] = ETLook.radiation.net_radiation_soil(od["rn_24"], od["sf_soil"])

    # del od["lai"]
    
    # **ETLook.resistance.soil_resistance***********************************************************

    if et_look_version == "dev":
        od["r_soil"] = ETLook.resistance.soil_resistance(od["se_root"], id["land_mask"], c.r_soil_pow, c.r_soil_min) #se_root was se_top
    elif et_look_version == "v2":
        od["r_soil"] = ETLook.resistance.soil_resistance(c.se_top, id["land_mask"], c.r_soil_pow, c.r_soil_min) #se_root was se_top

    # del od["se_root"]

    # **ETLook.resistance.soil_resistance***********************************************************

    od["ra_soil_init"] = ETLook.neutral.initial_soil_aerodynamic_resistance(id["u_24"], c.z_obs)

    # **ETLook.unstable.initial_friction_velocity_soil_daily***********************************************************

    od["u_star_24_soil_init"] = ETLook.unstable.initial_friction_velocity_soil_daily(od["u_b_24"], od["disp"], c.z_b)

    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************

    stc = ETLook.radiation.soil_thermal_conductivity(c.se_top)
    vhc = ETLook.radiation.volumetric_heat_capacity(c.se_top, c.porosity)
    dd = ETLook.radiation.damping_depth(stc, vhc)
    od["g0_bs"] = ETLook.radiation.bare_soil_heat_flux(doy, dd, stc, id["t_amp_year"], od["lat"])
    od["g0_24"] = ETLook.radiation.soil_heat_flux(od["g0_bs"], od["sf_soil"], id["land_mask"], od["rn_24_soil"], id["trans_24"], od["ra_24"], od["l_net"], c.rn_slope, c.rn_offset)
    od["e_24_init"] = ETLook.neutral.initial_daily_evaporation(od["rn_24_soil"], od["g0_24"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_soil"], od["ra_soil_init"])
    od["h_soil_24_init"] = ETLook.unstable.initial_sensible_heat_flux_soil_daily(od["rn_24_soil"], od["e_24_init"], od["g0_24"])

    # del od["sf_soil"], od["lat"], od["ra_soil_init"], od["g0_bs"], od["e_24_init"]
    
    # **ETLook.unstable.evaporation***********************************************************

    od["e_24"] = ETLook.unstable.evaporation(od["rn_24_soil"], od["g0_24"], od["ssvp_24"], od["ad_24"], od["vpd_24"], od["psy_24"], od["r_soil"], od["h_soil_24_init"], od["t_air_k_24"], od["u_star_24_soil_init"], od["disp"], od["u_b_24"], c.z_b, c.z_obs, c.iter_h)
    od["e_24_mm"] = ETLook.unstable.evaporation_mm(od["e_24"], od["lh_24"])
    od["et_24_mm"] = ETLook.evapotranspiration.et_actual_mm(od["e_24_mm"], od["t_24_mm"])
    
    if et_look_version == "dev":
        od["e_24_mm"][od["e_24_mm"] < 0] = 0
        od["et_24_mm"][od["et_24_mm"] < 0] = 0

    # del od["t_air_k_24"], od["u_b_24"], od["disp"], od["rn_24_soil"], od["r_soil"], od["u_star_24_soil_init"], od["h_soil_24_init"], od["e_24"], od["e_24_mm"]

    # **ETLook.unstable.evaporation***********************************************************
    od["rn_24_grass"] = ETLook.radiation.net_radiation_grass(od["ra_24"], od["l_net"], c.r0_grass)
    od["et_ref_24"] = ETLook.evapotranspiration.et_reference(od["rn_24_grass"], od["ad_24"], od["psy_24"], od["vpd_24"], od["ssvp_24"], id["u_24"])
    od["et_ref_24_mm"] = ETLook.evapotranspiration.et_reference_mm(od["et_ref_24"], od["lh_24"])

    if et_look_version == "v2":
        od["et_ref_24_mm"][od["et_ref_24_mm"] < 0] = 0

    # del od["vpd_24"], od["l_net"], od["ad_24"], od["psy_24"], od["ssvp_24"], od["rn_24_grass"], od["et_ref_24"], od["et_ref_24_mm"]
    
    if et_look_version == "dev":
        #ef_24 = ETLook.evapotranspiration.evaporative_fraction(et_24_mm, lh_24, rn_24, g0_24)    
        #eps_w = ETLook.stress.epsilon_soil_moisture(ef_24)    
        eps_a = ETLook.stress.epsilon_autotrophic_respiration()     
        od["lue"] = ETLook.biomass.lue(id["lue_max"], od["stress_temp"], od["stress_moist"], eps_a)
        od["fpar"] = ETLook.leaf.fpar(od["vc"], id["ndvi"])
        od["apar"] = ETLook.leaf.apar(od["ra_24"], od["fpar"])       
        od["biomass_prod"] = ETLook.biomass.biomass(od["apar"], od["lue"])     

    for var in vars.outputs.keys():
        if var in od.keys():
            od[var][np.isnan(QC)] = np.nan
            if vars.outputs[var]["output"]:
                # print("saving '{0}'.".format(vars.outputs[var]["file_path"]))
                PF.Save_as_tiff(vars.outputs[var]["file_path"], od[var], geo_ex, proj_ex)
        # else:
        #     print("'{0}' not found.".format(var))

    return
