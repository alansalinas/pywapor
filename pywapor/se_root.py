import xarray as xr
import os
import numpy as np
import pandas as pd
import pywapor
import pywapor.general as g
import pywapor.general.processing_functions as PF
import pywapor.et_look_dev as ETLook_dev
import pywapor.et_look_v2 as ETLook_v2
from pywapor.general.logger import log
from pywapor.general.compositer import calculate_ds
import copy

def main(input_data, se_root_version = "v2", export_vars = "default", export_to_tif = False):

    log.info("> SE_ROOT").add()

    # Version
    if se_root_version == "v2":
        ETLook = ETLook_v2
        log.info("--> Running SEroot_v2.")
    elif se_root_version == "dev":
        ETLook = ETLook_dev
        log.info("--> Running SEroot_dev.")
    c = ETLook.constants

    # Inputs
    if isinstance(input_data, str):
        ds = xr.open_dataset(input_data)
    else:
        ds = input_data
        input_data = ds.encoding["source"]

    if se_root_version == "dev":
        c.z0m_full = 0.04 + 0.01 * (ds.pixel_size - 30)/(250-30)
        ds["lst_zone_mean"] = lst_zone_mean(ds)

    ds["u_i"] = np.sqrt(ds["v2m_i"]**2 + ds["u2m_i"]**2)

    doy = [int(pd.Timestamp(x).strftime("%j")) for x in ds["time"].values]
    ds["doy"] = xr.DataArray(doy, coords = ds["time"].coords)
    dtime = [pd.Timestamp(x).hour + (pd.Timestamp(x).minute / 60) for x in ds["time"].values]
    ds["dtime"] = xr.DataArray(dtime, coords = ds["time"].coords)
    ds["sc"] = ETLook.solar_radiation.seasonal_correction(ds["doy"])
    ds["decl"] = ETLook.solar_radiation.declination(ds["doy"])
    ds["day_angle"] = ETLook.clear_sky_radiation.day_angle(ds["doy"])

    if "r0_bare" not in list(ds.variables): # TODO add this to pre_se_root and remove from pre_et_look
        ds["r0_bare"] = xr.ones_like(ds["ndvi"]) * c.r0_bare
    if "r0_full" not in list(ds.variables): # TODO add this to pre_se_root and remove from pre_et_look
        ds["r0_full"] = xr.ones_like(ds["ndvi"]) * c.r0_full

    if se_root_version == "dev":
        ds["t_air_i"] = xr.where(ds["t_air_i"] < -270, np.nan, ds["t_air_i"])

    ds["p_air_i"] = ETLook.meteo.air_pressure_kpa2mbar(ds["p_air_i"])
    ds["p_air_0_i"] = ETLook.meteo.air_pressure_kpa2mbar(ds["p_air_0_i"])

    ds["vc"] = ETLook.leaf.vegetation_cover(ds["ndvi"], c.nd_min, c.nd_max, c.vc_pow)

    ds["t_air_k_i"] = ETLook.meteo.air_temperature_kelvin_inst(ds["t_air_i"])
    ds["vp_i"] = ETLook.meteo.vapour_pressure_from_specific_humidity_inst(ds["qv_i"], ds["p_air_i"])
    ds["ad_moist_i"] = ETLook.meteo.moist_air_density_inst(ds["vp_i"], ds["t_air_k_i"])
    ds["ad_dry_i"] = ETLook.meteo.dry_air_density_inst(ds["p_air_i"], ds["vp_i"], ds["t_air_k_i"])
    ds["ad_i"] = ETLook.meteo.air_density_inst(ds["ad_dry_i"], ds["ad_moist_i"])
    ds["u_b_i_bare"] = ETLook.soil_moisture.wind_speed_blending_height_bare(ds["u_i"], c.z0m_bare, c.z_obs, c.z_b)
    ds["lon_rad"] = ETLook.solar_radiation.longitude_rad(ds["lon"])
    ds["lat_rad"] = ETLook.solar_radiation.latitude_rad(ds["lat"])
    ds["ha"] = ETLook.solar_radiation.hour_angle(ds["sc"], ds["dtime"], ds["lon_rad"])
    
    I0 = ETLook.clear_sky_radiation.solar_constant() # TODO this should come from ETLook.constants
    ds["ied"] = ETLook.clear_sky_radiation.inverse_earth_sun_distance(ds["day_angle"])
    ds["h0"] = ETLook.clear_sky_radiation.solar_elevation_angle(ds["lat_rad"], ds["decl"], ds["ha"])
    ds["h0ref"] = ETLook.clear_sky_radiation.solar_elevation_angle_refracted(ds["h0"])
    ds["m"] = ETLook.clear_sky_radiation.relative_optical_airmass(ds["p_air_i"], ds["p_air_0_i"], ds["h0ref"])
    ds["rotm"] = ETLook.clear_sky_radiation.rayleigh_optical_thickness(ds["m"])
    ds["Tl2"] = ETLook.clear_sky_radiation.linke_turbidity(ds["wv_i"], c.aod550_i, ds["p_air_i"], ds["p_air_0_i"])
    ds["G0"] = ETLook.clear_sky_radiation.extraterrestrial_irradiance_normal(I0, ds["ied"])
    ds["B0c"] = ETLook.clear_sky_radiation.beam_irradiance_normal_clear(ds["G0"], ds["Tl2"], ds["m"], ds["rotm"], ds["h0"])
    ds["Bhc"] = ETLook.clear_sky_radiation.beam_irradiance_horizontal_clear(ds["B0c"], ds["h0"])
    ds["Dhc"] = ETLook.clear_sky_radiation.diffuse_irradiance_horizontal_clear(ds["G0"], ds["Tl2"], ds["h0"])
    
    ds["ra_hor_clear_i"] = ETLook.clear_sky_radiation.ra_clear_horizontal(ds["Bhc"], ds["Dhc"])
    ds["emiss_atm_i"] = ETLook.soil_moisture.atmospheric_emissivity_inst(ds["vp_i"], ds["t_air_k_i"])
    
    ds["rn_bare"] = ETLook.soil_moisture.net_radiation_bare(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["lst"], ds["r0_bare"])
    ds["rn_full"] = ETLook.soil_moisture.net_radiation_full(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["lst"], ds["r0_full"])
    ds["h_bare"] = ETLook.soil_moisture.sensible_heat_flux_bare(ds["rn_bare"], c.fraction_h_bare)
    ds["h_full"] = ETLook.soil_moisture.sensible_heat_flux_full(ds["rn_full"], c.fraction_h_full)
    ds["u_b_i_full"] = ETLook.soil_moisture.wind_speed_blending_height_full_inst(ds["u_i"], c.z0m_full, c.z_obs, c.z_b)
    
    ds["u_star_i_bare"] = ETLook.soil_moisture.friction_velocity_bare_inst(ds["u_b_i_bare"], c.z0m_bare, c.disp_bare, c.z_b)
    ds["u_star_i_full"] = ETLook.soil_moisture.friction_velocity_full_inst(ds["u_b_i_full"], c.z0m_full, c.disp_full, c.z_b)
    ds["L_bare"] = ETLook.soil_moisture.monin_obukhov_length_bare(ds["h_bare"], ds["ad_i"], ds["u_star_i_bare"], ds["t_air_k_i"])
    ds["L_full"] = ETLook.soil_moisture.monin_obukhov_length_full(ds["h_full"], ds["ad_i"], ds["u_star_i_full"], ds["t_air_k_i"])
    
    ds["u_i_soil"] = ETLook.soil_moisture.wind_speed_soil_inst(ds["u_i"], ds["L_bare"], c.z_obs)
    ds["ras"] = ETLook.soil_moisture.aerodynamical_resistance_soil(ds["u_i_soil"])
    ds["raa"] = ETLook.soil_moisture.aerodynamical_resistance_bare(ds["u_i"], ds["L_bare"], c.z0m_bare, c.disp_bare, c.z_obs)
    ds["rac"] = ETLook.soil_moisture.aerodynamical_resistance_full(ds["u_i"], ds["L_full"], c.z0m_full, c.disp_full, c.z_obs)
    
    ds["t_max_bare"] = ETLook.soil_moisture.maximum_temperature_bare(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["raa"], ds["ras"], ds["r0_bare"])
    ds["t_max_full"] = ETLook.soil_moisture.maximum_temperature_full(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["rac"], ds["r0_full"])

    # od["w_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])
    # od["t_dew_i"] = ETLook.soil_moisture.dew_point_temperature_inst(od["vp_i"])
    # od["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst(id["t_air_i"], od["t_dew_i"], id["p_air_i"])
    ds["t_wet_i"] = ETLook.soil_moisture.wet_bulb_temperature_inst_new(ds["t_air_i"], ds["qv_i"], ds["p_air_i"])
    #t_wet_i_new2 = ETLook.soil_moisture.wet_bulb_temperature_inst_new2(t_air_i)

    ds["lst_max"] = ETLook.soil_moisture.maximum_temperature(ds["t_max_bare"], ds["t_max_full"], ds["vc"])

    if se_root_version == "v2":
        ds["t_wet_k_i"] = ETLook.meteo.wet_bulb_temperature_kelvin_inst(ds["t_wet_i"])
        ds["lst_min"] = ETLook.soil_moisture.minimum_temperature(ds["t_wet_k_i"], ds["t_air_k_i"], ds["vc"])
    elif se_root_version == "dev":
        ds["t_min_bare"] = ETLook.soil_moisture.minimum_temperature_bare(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["raa"], ds["ras"], ds["lst_zone_mean"], c.r0_bare_wet)
        ds["t_min_full"] = ETLook.soil_moisture.minimum_temperature_full(ds["ra_hor_clear_i"], ds["emiss_atm_i"], ds["t_air_k_i"], ds["ad_i"], ds["rac"], ds["lst_zone_mean"], ds["r0_full"])  
        ds["lst_min"] = ETLook.soil_moisture.maximum_temperature(ds["t_min_bare"], ds["t_min_full"], ds["vc"])

    ds["se_root"] = ETLook.soil_moisture.soil_moisture_from_maximum_temperature(ds["lst_max"], ds["lst"], ds["lst_min"])

    if export_vars == "all":
        ...
    elif export_vars == "default":
        keep_vars = ['se_root']
        ds = PF.ds_remove_except(ds, keep_vars)
    elif isinstance(export_vars, list):
        keep_vars = copy.copy(export_vars)
        ds = PF.ds_remove_except(ds, keep_vars)
    else:
        raise ValueError

    ds = g.variables.fill_attrs(ds)

    fp, fn = os.path.split(input_data)
    fn = fn.replace("_input", "_output")
    ds, fh = calculate_ds(ds, os.path.join(fp, fn), "--> Saving outputs.")

    log.sub().info("< SE_ROOT")

    if export_to_tif:
        files = PF.export_ds_to_tif(ds, keep_vars, None)
        ds.close()
        # os.remove(fh)
        return files
    else:
        return ds

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

def lst_zone_mean(ds):
    # TODO This function needs to be replaced by something like pywapor.enhancers.temperature.local_mean

    geo_ex = ds.geotransform
    proj_ex = ds.projection

    size_y, size_x = ds["ndvi"].isel(time = 0).shape
    size_y_zone = int(np.ceil(size_y/200))
    size_x_zone = int(np.ceil(size_x/200)) 
    array_fake = np.ones([size_y_zone, size_x_zone])
    geo_new = tuple([geo_ex[0], geo_ex[1] * 200, geo_ex[2], geo_ex[3], geo_ex[4], geo_ex[5]*200])
    MEM_file = PF.Save_as_MEM(array_fake, geo_new, proj_ex)

    lst_zone_mean_full = np.ones_like(ds["ndvi"] * np.nan)

    for i, t in enumerate(ds.time):

        lst = ds["lst"].sel(time = t).values
        lst_filename = PF.Save_as_MEM(lst, geo_new, proj_ex)
        
        dest_lst_zone_large = PF.reproject_dataset_example(lst_filename, MEM_file, 4)
        lst_zone_mean_large = dest_lst_zone_large.GetRasterBand(1).ReadAsArray()
        lst_zone_mean_large[lst_zone_mean_large==0] = -9999
        lst_zone_mean_large[np.isnan(lst_zone_mean_large)] = -9999

        if np.nanmax(lst_zone_mean_large) == -9999:
            for x in range(0, size_x_zone):
                for y in range(0, size_y_zone):
                    lst_zone_mean_large[y, x] = np.nanmean(lst[y*200:np.minimum((y+1)*200, size_y-1), x*200:np.minimum((x+1)*200, size_x-1)])
            lst_zone_mean_large[np.isnan(lst_zone_mean_large)] = -9999

        lst_zone_mean_large = PF.gap_filling(lst_zone_mean_large, -9999, 1)
        dest_lst_zone_large = PF.Save_as_MEM(lst_zone_mean_large, geo_new, proj_ex)
        
        dest_lst_zone = PF.reproject_dataset_example(dest_lst_zone_large, lst_filename, 6)
        lst_zone_mean = dest_lst_zone.GetRasterBand(1).ReadAsArray()

        lst_zone_mean_full[i, ...] = lst_zone_mean
        
    lst_zone_mean_da = xr.DataArray(lst_zone_mean_full, coords = ds.ndvi.coords)

    return lst_zone_mean_da

if __name__ == "__main__":

    import pywapor

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]

    startdate = "2021-07-06"
    enddate = "2021-07-07"
    composite_length = 1
    level = "level_2"
    extra_sources = None
    extra_source_locations = None

    my_custom_level = {
        # Main inputs
        "ndvi":         ["PROBAV", "MOD13", "MYD13"],
        "r0":           ["PROBAV", "MCD43"],
        "lst":          ["MOD11", "MYD11"],
        "lulc":         ["GLOBCOVER"],
        "z":            ["SRTM"],
        "p_24":         ["CHIRPS"],
        "ra_24":        ["MERRA2"],

        # Daily meteo 
        't_air_24':     ["MERRA2", "GEOS5"],
        't_air_min_24': ["MERRA2"], 
        't_air_max_24': ["MERRA2"],
        'u2m_24':       ["GEOS5"],
        'v2m_24':       ["GEOS5"],
        'p_air_0_24':   ["MERRA2"],
        'qv_24':        ["MERRA2", "GEOS5"],

        # Instanteneous meteo
        "t_air_i":      ["MERRA2"],
        "u2m_i":        ["MERRA2"],
        "v2m_i":        ["MERRA2"],
        "qv_i":         ["MERRA2"],
        "wv_i":         ["MERRA2"],
        "p_air_i":      ["MERRA2"],
        "p_air_0_i":    ["MERRA2"],

        # Temporal constants
        "lw_offset":    ["STATICS"],
        "lw_slope":     ["STATICS"],
        "r0_bare":      ["STATICS"],
        "r0_full":      ["STATICS"],
        "rn_offset":    ["STATICS"],
        "rn_slope":     ["STATICS"],
        "t_amp_year":   ["STATICS"],
        "t_opt":        ["STATICS"],
        "vpd_slope":    ["STATICS"],
        "z_oro":        ["STATICS"],

        # Level name
        "name": "test",
    }

    diagnostics = None

    ds, fh = pywapor.pre_se_root.main(project_folder, startdate, enddate, latlim, 
                                        lonlim, level = my_custom_level,
                                        extra_sources = extra_sources, 
                                        extra_source_locations = extra_source_locations)

    se_root_version = "v2"

    ds_out = main(ds, se_root_version = se_root_version, export_vars = "all", export_to_tif = False)


    