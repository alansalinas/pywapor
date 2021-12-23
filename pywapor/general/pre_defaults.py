
def composite_defaults():

    cdefaults = {

        'NDVI': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "ndvi",
                    "var_unit": "-",
                },
        'ALBEDO': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "r0",
                    "var_unit": "-",
                },
        'LULC': {
                    "composite_type": 'mean',
                    "temporal_interp": 'nearest',
                    "temporal_interp_freq": 1,
                    "spatial_interp": "linear",
                    "var_name": "lulc",
                    "var_unit": "-",
                },
        'DEM': {
                    "composite_type": False,
                    "temporal_interp": False,
                    "spatial_interp": "linear",
                    "temporal_interp_freq": 1,
                    "var_name": "z",
                    "var_unit": "m",
                },
        'PRECIPITATION': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "P_24",
                    "var_unit": "mm/day",
                },
        'SE_ROOT': {
                    "composite_type": 0.85,
                    "temporal_interp": False,
                    "temporal_interp_freq": "2D",
                    "spatial_interp": "nearest",
                    "var_name": "se_root",
                    "var_unit": "-",
                },
        'SOLAR_RADIATION': {
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
        "u2m_24":"m/s", "v2m_24":"m/s","qv_24":"kg/kg","p_air_24_0":"kPa"}

    for name, unit in meteo_names_units.items():
        cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **meteo}

    return cdefaults