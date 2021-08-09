inputs = {
    "ALBEDO":                   {"time": "daily", "array_name": "r0"},
    "NDVI":                     {"time": "daily", "array_name": "ndvi"}, 
    "LST":                      {"time": "daily", "array_name": "lst"},
    "Time":                     {"time": "daily", "array_name": "dtime"},
    "Lat":                      {"time": "static", "array_name": "lat_deg"},
    "Lon":                      {"time": "static", "array_name": "lon_deg"},
    "DEM":                      {"time": "static", "array_name": "z"},
    "Slope":                    {"time": "static", "array_name": "slope_deg"},
    "Aspect":                   {"time": "static", "array_name": "aspect_deg"},
    "LandMask":                 {"time": "yearly", "array_name": "land_mask"},
    "LUEmax":                   {"time": "yearly", "array_name": "lue_max"},
    "Pair_24_0":                {"time": "daily", "array_name": "p_air_0_24"},
    "Pair_inst_0":              {"time": "daily", "array_name": "p_air_0_i"},
    "Pair_inst":                {"time": "daily", "array_name": "p_air_i"},
    "Precipitation":            {"time": "daily", "array_name": "P_24"},
    "qv_24":                    {"time": "daily", "array_name": "qv_24"},
    "qv_inst":                  {"time": "daily", "array_name": "qv_i"},
    "tair_24":                  {"time": "daily", "array_name": "t_air_24"},
    "tair_max_24":              {"time": "daily", "array_name": "t_air_max_24"},
    "tair_min_24":              {"time": "daily", "array_name": "t_air_min_24"},
    "tair_inst":                {"time": "daily", "array_name": "t_air_i"},
    "Tair_amp":                 {"time": "yearly", "array_name": "t_amp_year"},
    "wind_24":                  {"time": "daily", "array_name": "u_24"},
    "wind_inst":                {"time": "daily", "array_name": "u_i"},
    "wv_inst":                  {"time": "daily", "array_name": "wv_i"},
    "Trans_24":                 {"time": "daily", "array_name": "trans_24"},
    "se_root":                  {"time": "daily", "array_name": "se_root"},
    "Maximum_Obstacle_Height":  {"time": "yearly", "array_name": "z_obst_max"},
    "Bulk_Stomatal_resistance": {"time": "yearly", "array_name": "rs_min"}
}

outputs = {
    'B0c':                  {'file_name': 'B0c'},
    'Bhc':                  {'file_name': 'Bhc'},
    'Dhc':                  {'file_name': 'Dhc'},
    'Tl2':                  {'file_name': 'Tl2'},
    'ad_24':                {'file_name': 'ad_24'},
    'ad_dry_24':            {'file_name': 'ad_dry_24'},
    'ad_dry_i':             {'file_name': 'ad_dry_i'},
    'ad_i':                 {'file_name': 'ad_i'},
    'ad_moist_24':          {'file_name': 'ad_moist_24'},
    'ad_moist_i':           {'file_name': 'ad_moist_i'},
    'aspect':               {'file_name': 'aspect'},
    'biomass_prod':         {'file_name': 'biomass_prod_kg-ha'}, #
    'diffusion_index':      {'file_name': 'diffusion_index'},
    'disp':                 {'file_name': 'disp'},
    'e_24':                 {'file_name': 'e_24'},
    'e_24_init':            {'file_name': 'e_24_init'},
    'e_24_mm':              {'file_name': 'e_24_mm'},
    'emiss_atm_i':          {'file_name': 'emiss_atm_i'},
    'et_24_mm':             {'file_name': 'et_24_mm'},
    'et_ref_24':            {'file_name': 'et_ref_24'},
    'et_ref_24_mm':         {'file_name': 'et_ref_24_mm'},
    'g0_24':                {'file_name': 'g0_24'},
    'g0_bs':                {'file_name': 'g0_bs'},
    'h0':                   {'file_name': 'h0'},
    'h0ref':                {'file_name': 'h0ref'},
    'h_canopy_24_init':     {'file_name': 'h_canopy_24_init'},
    'h_soil_24_init':       {'file_name': 'h_soil_24_init'},
    'ha':                   {'file_name': 'ha'},
    'int_mm':               {'file_name': 'int_mm'},
    'int_wm2':              {'file_name': 'int_wm2'},
    'l_net':                {'file_name': 'l_net'},
    'lai':                  {'file_name': 'LAI'}, #
    'lai_eff':              {'file_name': 'LAI_eff'}, #
    'lat':                  {'file_name': 'lat'},
    'lh_24':                {'file_name': 'lh_24'},
    'lon':                  {'file_name': 'lon'},
    'lst_max':              {'file_name': 'lst_max'},
    'lst_min':              {'file_name': 'lst_min'},
    'm':                    {'file_name': 'm'},
    'p_air_24':             {'file_name': 'p_air_24'},
    'psy_24':               {'file_name': 'psy_24'},
    'r_canopy':             {'file_name': 'r_canopy'},
    'r_canopy_0':           {'file_name': 'r_canopy_0'},
    'r_soil':               {'file_name': 'r_soil'},
    'ra_24':                {'file_name': 'ra_24'},
    'ra_24_toa':            {'file_name': 'ra_24_toa'},
    'ra_24_toa_flat':       {'file_name': 'ra_24_toa_flat'}, #
    'ra_canopy_init':       {'file_name': 'ra_canopy_init'},
    'ra_hor_clear_i':       {'file_name': 'ra_hor_clear_i'},
    'ra_soil_init':         {'file_name': 'ra_soil_init'},
    'raa':                  {'file_name': 'raa'},
    'rac':                  {'file_name': 'rac'},
    'ras':                  {'file_name': 'ras'},
    'rn_24':                {'file_name': 'rn_24'},
    'rn_24_canopy':         {'file_name': 'rn_24_canopy'},
    'rn_24_grass':          {'file_name': 'rn_24_grass'},
    'rn_24_soil':           {'file_name': 'rn_24_soil'},
    'rn_bare':              {'file_name': 'rn_bare'},
    'rn_full':              {'file_name': 'rn_full'},
    'rotm':                 {'file_name': 'rotm'},
    'se_root':              {'file_name': 'se_root'},
    'sf_soil':              {'file_name': 'sf_soil'},
    'slope':                {'file_name': 'slope'},
    'ssvp_24':              {'file_name': 'ssvp_24'},
    'stress_moist':         {'file_name': 'stress_moist'},
    'stress_rad':           {'file_name': 'stress_rad'},
    'stress_temp':          {'file_name': 'stress_temp'},
    'stress_vpd':           {'file_name': 'stress_vpd'},
    'svp_24':               {'file_name': 'svp_24'},
    't_24':                 {'file_name': 't_24'},
    't_24_init':            {'file_name': 't_24_init'},
    't_24_mm':              {'file_name': 't_24_mm'},
    'tpot_24':              {'file_name': 'tpot_24'},
    'tpot_24_mm':           {'file_name': 'tpot_24_mm'},
    't_air_k_24':           {'file_name': 't_air_k_24'},
    't_air_k_i':            {'file_name': 't_air_k_i'},
    't_dew_i':              {'file_name': 't_dew_i'},
    't_max_bare':           {'file_name': 't_max_bare'},
    't_max_full':           {'file_name': 't_max_full'},
    't_min_bare':           {'file_name': 't_min_bare'},
    't_min_full':           {'file_name': 't_min_full'},
    't_wet_i':              {'file_name': 't_wet_i'},
    't_wet_k_i':            {'file_name': 't_wet_k_i'},
    'u_b_24':               {'file_name': 'u_b_24'},
    'u_b_i_bare':           {'file_name': 'u_b_i_bare'},
    'u_b_i_full':           {'file_name': 'u_b_i_full'},
    'u_i_soil':             {'file_name': 'u_i_soil'},
    'u_star_24_init':       {'file_name': 'u_star_24_init'},
    'u_star_24_soil_init':  {'file_name': 'u_star_24_soil_init'},
    'u_star_i_bare':        {'file_name': 'u_star_i_bare'},
    'u_star_i_full':        {'file_name': 'u_star_i_full'},
    'vc':                   {'file_name': 'vc'},
    'vp_24':                {'file_name': 'vp_24'},
    'vp_i':                 {'file_name': 'vp_i'},
    'vpd_24':               {'file_name': 'vpd_24'},
    'w_i':                  {'file_name': 'w_i'},
    'ws':                   {'file_name': 'ws'},
    'z0m':                  {'file_name': 'z0m'},
    'z_obst':               {'file_name': 'z_obst'},
    'z_oro':                {'file_name': 'z_oro'},
    'L_bare':               {'file_name': 'L_bare'},
    'L_full':               {'file_name': 'L_full'},
}

def get_temp_input_data_reqs():
    temporal_input_data_req = [
        {"name": "Albedo",
        "unit": "[-]",
        "quantity": "Albedo",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/ALBEDO_{year}{month}{day}.tif"},

        {"name": "Land Surface Temperature",
        "unit": "[K]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/LST_{year}{month}{day}.tif"},

        {"name": "Normalized Difference Vegetation Index",
        "unit": "[-]",
        "quantity": "Normalized Difference Vegetation Index",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/NDVI_{year}{month}{day}.tif"},

        {"name": "Air Pressure at sea level (daily average)",
        "unit": "[kPa]",
        "quantity": "Air Pressure",
        "level": "sea",
        "time": "daily average",
        "filepath": "{input_folder}/{year}{month}{day}/Pair_24_0_{year}{month}{day}.tif"},

        {"name": "Air Pressure at sea level (instanteneous)", 
        "unit": "[kPa]",
        "quantity": "Air Pressure",
        "level": "sea",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/Pair_inst_0_{year}{month}{day}.tif"},

        {"name": "Air Pressure at surface level (instanteneous)",
        "unit": "[kPa]",
        "quantity": "Air Pressure",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/Pair_inst_{year}{month}{day}.tif"},

        {"name": "Precipitation",
        "unit": "[mm/day]",
        "quantity": "Precipitation",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/Precipitation_{year}{month}{day}.tif"},

        {"name": "Specific Humidity (daily average)", 
        "unit": "[kg/kg]",
        "quantity": "Specific Humidity",
        "level": "surface",
        "time": "daily average",
        "filepath": "{input_folder}/{year}{month}{day}/qv_24_{year}{month}{day}.tif"},

        {"name": "Specific Humidity (instanteneous)", "unit": "[kg/kg]",
        "quantity": "Specific Humidity",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/qv_inst_{year}{month}{day}.tif"},

        {"name": "Air Temperature (daily average)", "unit": "[C]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "daily average",
        "filepath": "{input_folder}/{year}{month}{day}/tair_24_{year}{month}{day}.tif"},

        {"name": "Air Temperature (instanteneous)", "unit": "[C]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/tair_inst_{year}{month}{day}.tif"},

        {"name": "Air Temperature (daily maximum)", "unit": "[C]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "daily maximum",
        "filepath": "{input_folder}/{year}{month}{day}/tair_max_24_{year}{month}{day}.tif"},

        {"name": "Air Temperature (daily minimum)", "unit": "[C]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "daily minimum",
        "filepath": "{input_folder}/{year}{month}{day}/tair_min_24_{year}{month}{day}.tif"},

        {"name": "Transmissivity", "unit": "[-]",
        "quantity": "Transmissivity",
        "level": "column",
        "time": "daily average",
        "filepath": "{input_folder}/{year}{month}{day}/Trans_24_{year}{month}{day}.tif"},

        {"name": "Windspeed (daily average)", "unit": "[m/s]",
        "quantity": "Windspeed",
        "level": "surface",
        "time": "daily average",
        "filepath": "{input_folder}/{year}{month}{day}/wind_24_{year}{month}{day}.tif"},

        {"name": "Windspeed (instanteneous)", "unit": "[m/s]",
        "quantity": "Windspeed",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/wind_inst_{year}{month}{day}.tif"},

        {"name": "Total Precipitable Water Vapour", "unit": "[mm]",
        "quantity": "Total Precipitable Water Vapour",
        "level": "column",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/wv_inst_{year}{month}{day}.tif"},

        {"name": "Instantaneous Data Time", "unit": "[hour]",
        "quantity": "Time",
        "level": "surface",
        "time": "instanteneous",
        "filepath": "{input_folder}/{year}{month}{day}/Time_{year}{month}{day}.tif"},
    ]

    return temporal_input_data_req

def get_static_input_data_reqs():
    static_input_data_req = [
        {"name": "Latitude",
        "unit": "[DD]",
        "quantity": "Latitude",
        "level": "surface",
        "time": "invariant",
        "filepath": "{input_folder}/Static/Lat.tif"},

        {"name": "Longitude",
        "unit": "[DD]",
        "quantity": "Longitude",
        "level": "surface",
        "time": "invariant",
        "filepath": "{input_folder}/Static/Lon.tif"},

        {"name": "Slope",
        "unit": "[degrees]",
        "quantity": "Slope",
        "level": "surface",
        "time": "invariant",
        "filepath": "{input_folder}/Static/Slope.tif"},

        {"name": "Slope Aspect",
        "unit": "[degrees]",
        "quantity": "Aspect",
        "level": "surface",
        "time": "invariant",
        "filepath": "{input_folder}/Static/Aspect.tif"},

        {"name": "Bulk Stomatal Resistance",
        "unit": "[s/m]",
        "quantity": "Resistance",
        "level": "surface",
        "time": "yearly",
        "filepath": "{input_folder}/Static/Bulk_Stomatal_resistance_{year}.tif"},

        {"name": "Digital Elevation Model",
        "unit": "[m.a.s.l]",
        "quantity": "Altitude",
        "level": "surface",
        "time": "invariant",
        "filepath": "{input_folder}/Static/DEM.tif"},

        {"name": "Landmask",
        "unit": "[-]",
        "quantity": "Landmask",
        "level": "surface",
        "time": "yearly",
        "filepath": "{input_folder}/Static/LandMask_{year}.tif"},

        {"name": "Maximum Light Use Efficiency",
        "unit": "[gr/MJ]",
        "quantity": "Efficiency",
        "level": "surface",
        "time": "yearly",
        "filepath": "{input_folder}/Static/LUEmax_{year}.tif"},

        {"name": "Maximum Obstacle Height",
        "unit": "[m]",
        "quantity": "Height",
        "level": "surface",
        "time": "yearly",
        "filepath": "{input_folder}/Static/Maximum_Obstacle_Height_{year}.tif"},

        {"name": "Air Temperature (yearly amplitude)",
        "unit": "[K]",
        "quantity": "Temperature",
        "level": "surface",
        "time": "yearly",
        "filepath": "{input_folder}/Static/Tair_amp_{year}.tif"},

    ]
    return static_input_data_req