from pywapor.enhancers.temperature import lapse_rate
from functools import partial
import pywapor.enhancers.dem as dem
from pywapor.enhancers.temperature import kelvin_to_celsius
import pywapor.enhancers.lulc as lulc
from functools import partial

def composite_enhancements_defaults():
    composite_enhancements = {
        "t_air_24":     [lapse_rate],
        "t_air_min_24": [lapse_rate],
        "t_air_max_24": [lapse_rate],
        "z":            [partial(dem.to_slope, out_var = "slope"), 
                        partial(dem.to_aspect, out_var = "aspect"),
                        partial(dem.to_lat, out_var = "lat_deg"),
                        partial(dem.to_lon, out_var = "lon_deg"),]
    }
    return composite_enhancements

def source_enhancements_defaults():

    def remove_var(ds, var):
        return ds.drop_vars([var])

    source_enhancements = {
        ("MERRA2",  "t_air_24"):        [kelvin_to_celsius],
        ("MERRA2",  "t_air_min_24"):    [kelvin_to_celsius],
        ("MERRA2",  "t_air_max_24"):    [kelvin_to_celsius],
        ("GEOS5",   "t_air_24"):        [kelvin_to_celsius],
        ("GEOS5",   "t_air_min_24"):    [kelvin_to_celsius],
        ("GEOS5",   "t_air_max_24"):    [kelvin_to_celsius],
        ("GLOBCOVER", "lulc"):          [partial(lulc.lulc_to_x, out_var = "land_mask", 
                                                    convertor = lulc.globcover_to_land_mask()),
                                        partial(lulc.lulc_to_x, out_var = "rs_min", 
                                                    convertor = lulc.globcover_to_rs_min()),
                                        partial(lulc.lulc_to_x, out_var = "lue_max", 
                                                    convertor = lulc.globcover_to_lue_max()),
                                        partial(lulc.lulc_to_x, out_var = "z_obst_max", 
                                                    convertor = lulc.globcover_to_z_obst_max()),
                                        remove_var,
                                        ],
        ("WAPOR", "lulc"):              [partial(lulc.lulc_to_x, out_var = "land_mask", 
                                                    convertor = lulc.wapor_to_land_mask()),
                                        partial(lulc.lulc_to_x, out_var = "rs_min", 
                                                    convertor = lulc.wapor_to_rs_min()),
                                        partial(lulc.lulc_to_x, out_var = "lue_max", 
                                                    convertor = lulc.wapor_to_lue_max()),
                                        partial(lulc.lulc_to_x, out_var = "z_obst_max", 
                                                    convertor = lulc.wapor_to_z_obst_max()),
                                        remove_var,
                                        ],
    }
    return source_enhancements

def composite_defaults():

    cdefaults = {

        'ndvi': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "spatial_interp": "nearest",
                    "var_name": "ndvi",
                    "var_unit": "-",
                },
        'r0': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "spatial_interp": "nearest",
                    "var_name": "r0",
                    "var_unit": "-",
                },
        'lulc': {
                    "composite_type": 'mean',
                    "temporal_interp": 'nearest',
                    "spatial_interp": "linear",
                    "var_name": "lulc",
                    "var_unit": "-",
                },
        'z': {
                    "composite_type": False,
                    "temporal_interp": False,
                    "spatial_interp": "linear",
                    "var_name": "z",
                    "var_unit": "m",
                },
        'p_24': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "spatial_interp": "nearest",
                    "var_name": "p_24",
                    "var_unit": "mm/day",
                },
        'se_root': {
                    "composite_type": "max",
                    "temporal_interp": False,
                    "spatial_interp": "nearest",
                    "var_name": "se_root",
                    "var_unit": "-",
                },
        'ra_24': {
                    "composite_type": "mean",
                    "temporal_interp": "linear",
                    "spatial_interp": "linear",
                    "var_name": "ra_24",
                    "var_unit": "-",
                }
        }

    meteo = {
            "composite_type": "mean",
            "temporal_interp": "linear",
            "spatial_interp": "linear",
            }

    meteo_names_units = {"t_air_max_24": "°C", "t_air_min_24":"°C", "t_air_24":"°C", 
        "u2m_24":"m/s", "v2m_24":"m/s","qv_24":"kg/kg","p_air_0_24":"kPa"}

    for name, unit in meteo_names_units.items():
        cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **meteo}

    static = {
            "composite_type": False,
            "temporal_interp": False,
            "spatial_interp": "linear",    
    }

    statics_names_units = {'land_mask': "-",
                        'lw_offset': "-",
                        'lw_slope': "-",
                        'r0_bare': "-",
                        'r0_full': "-",
                        'rn_offset': "-",
                        'rn_slope': "-",
                        'rs_min': "sm-1",
                        't_amp_year': "C", 
                        't_opt': "C",
                        'vpd_slope': "mbar-1", 
                        'z_obst_max': "m",
                        'z_oro': "m",
    }

    for name, unit in statics_names_units.items():
        cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **static}

    return cdefaults