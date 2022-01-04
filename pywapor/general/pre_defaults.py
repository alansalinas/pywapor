
def composite_defaults():

    cdefaults = {

        'ndvi': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "ndvi",
                    "var_unit": "-",
                },
        'r0': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "r0",
                    "var_unit": "-",
                },
        'lulc': {
                    "composite_type": 'mean',
                    "temporal_interp": 'nearest',
                    "temporal_interp_freq": 1,
                    "spatial_interp": "linear",
                    "var_name": "lulc",
                    "var_unit": "-",
                },
        'z': {
                    "composite_type": False,
                    "temporal_interp": False,
                    "spatial_interp": "linear",
                    "temporal_interp_freq": 1,
                    "var_name": "z",
                    "var_unit": "m",
                },
        'p_24': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "p_24",
                    "var_unit": "mm/day",
                },
        'se_root': {
                    "composite_type": 0.85,
                    "temporal_interp": False,
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "se_root",
                    "var_unit": "-",
                },
        'ra_24': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "linear",
                    "var_name": "ra_24",
                    "var_unit": "-",
                }
        }

    meteo = {
            "composite_type": "mean",
            "temporal_interp": "linear",
            "temporal_interp_freq": "2D",
            "spatial_interp": "linear",
            }

    meteo_names_units = {"t_air_max_24": "°C", "t_air_min_24":"°C", "t_air_24":"°C", 
        "u2m_24":"m/s", "v2m_24":"m/s","qv_24":"kg/kg","p_air_0_24":"kPa"}

    for name, unit in meteo_names_units.items():
        cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **meteo}

    static = {
            "composite_type": False,
            "temporal_interp": False,
            "temporal_interp_freq": 1,
            "spatial_interp": "linear",    
    }

    statics_names_units = {'land_mask': "-",
                        'lw_offset': "-",
                        'lw_slope': "-",
                        'r0_bare': "-",
                        'r0_full': "-",
                        'rn_offset': "-",
                        'rn_slope': "-",
                        'rs_min': "s/m",
                        't_amp_year': "°C", 
                        't_opt': "°C",
                        'vpd_slope': "1/mbar", 
                        'z_obst_max': "m",
                        'z_oro': "m",
    }

    for name, unit in statics_names_units.items():
        cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **static}

    return cdefaults