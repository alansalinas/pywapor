"""Module containing functions that generate some default settings for `pywapor.pre_et_look`
and `pywapor.pre_se_root`.
"""
from pywapor.enhancers.temperature import lapse_rate
from functools import partial
import pywapor.enhancers.dem as dem
from pywapor.enhancers.temperature import kelvin_to_celsius
import pywapor.enhancers.lulc as lulc
from functools import partial

def constants_defaults():
    """Return a dictionary with default parameters.

    Returns
    -------
    dict
        Default constants.
    """

    c = {
        'nd_min': 0.125,
        'nd_max': 0.8,
        'vc_pow': 0.7,
        'vc_min': 0,
        'vc_max': 0.9677324224821418,
        'lai_pow': -0.45,
        'diffusion_slope': -1.33,
        'diffusion_intercept': 1.15,
        't_min': 0,
        't_max': 50,
        'vp_slope': 0.14,
        'vp_offset': 0.34,
        'int_max': 0.2,
        'tenacity': 1.5,
        'rcan_max': 1000000,
        'ndvi_obs_min': 0.25,
        'ndvi_obs_max': 0.75,
        'obs_fr': 0.25,
        'z_obs': 2,
        'z_b': 100,
        'c1': 1,
        'iter_h': 3,
        'r_soil_pow': -2.1,
        'r_soil_min': 800,
        'se_top': 0.5,
        'porosity': 0.4,
        'r0_grass': 0.23,
        'eps_a': 0.5,

        # se_root
        'z0m_full': 0.1,
        'z0m_bare': 0.001,
        'aod550_i': 0.01,
        'fraction_h_bare': 0.65,
        'fraction_h_full': 0.95,
        'disp_bare': 0.0,
        'disp_full': 0.667,
        'r0_bare_wet': 0.2,
        'IO': 1367.,

        # Biomass
        'dh_ap': 52750, 
        'd_s': 704.98, 
        'dh_dp': 211000,
        'ar_slo': 0.0, 
        'ar_int': 0.5,
        'fpar_slope': 1.257, 
        'fpar_offset': -0.161,
        'o2': 20.9,
        'co2_ref': 281,
        'gcgdm': 0.4,
        'phot_eff': 2.49,

    }

    return c


# def composite_enhancements_defaults():
#     """Returns a dictionary with the default functions to apply to the calculated
#     composites.

#     Returns
#     -------
#     dict
#         Keys are variable names, values are lists of functions.
#     """
#     composite_enhancements = {
#         "t_air_24":     [lapse_rate],
#         "t_air_min_24": [lapse_rate],
#         "t_air_max_24": [lapse_rate],
#         "z":            [partial(dem.to_slope, out_var = "slope"), 
#                         partial(dem.to_aspect, out_var = "aspect"),
#                         partial(dem.to_lat, out_var = "lat_deg"),
#                         partial(dem.to_lon, out_var = "lon_deg"),]
#     }
#     return composite_enhancements

# def source_enhancements_defaults():
#     """Returns a dictionary with the default functions to apply to source
#     data.

#     Returns
#     -------
#     dict
#         Keys are tuples like `("source", "variable")`, values are lists of functions.
#     """

#     def remove_var(ds, var):
#         return ds.drop_vars([var])

#     source_enhancements = {
#         ("MERRA2",  "t_air_24"):        [kelvin_to_celsius],#
#         ("MERRA2",  "t_air_min_24"):    [kelvin_to_celsius],#
#         ("MERRA2",  "t_air_max_24"):    [kelvin_to_celsius],#
#         ("GEOS5",   "t_air_24"):        [kelvin_to_celsius],#
#         ("GEOS5",   "t_air_min_24"):    [kelvin_to_celsius],#
#         ("GEOS5",   "t_air_max_24"):    [kelvin_to_celsius],#
#         ("GLOBCOVER", "lulc"):          [partial(lulc.lulc_to_x, out_var = "land_mask", 
#                                                     convertor = lulc.globcover_to_land_mask()),
#                                         partial(lulc.lulc_to_x, out_var = "rs_min", 
#                                                     convertor = lulc.globcover_to_rs_min()),
#                                         partial(lulc.lulc_to_x, out_var = "lue_max", 
#                                                     convertor = lulc.globcover_to_lue_max()),
#                                         partial(lulc.lulc_to_x, out_var = "z_obst_max", 
#                                                     convertor = lulc.globcover_to_z_obst_max()),
#                                         remove_var,
#                                         ],
#         ("WAPOR", "lulc"):              [partial(lulc.lulc_to_x, out_var = "land_mask", 
#                                                     convertor = lulc.wapor_to_land_mask()),
#                                         partial(lulc.lulc_to_x, out_var = "rs_min", 
#                                                     convertor = lulc.wapor_to_rs_min()),
#                                         partial(lulc.lulc_to_x, out_var = "lue_max", 
#                                                     convertor = lulc.wapor_to_lue_max()),
#                                         partial(lulc.lulc_to_x, out_var = "z_obst_max", 
#                                                     convertor = lulc.wapor_to_z_obst_max()),
#                                         remove_var,
#                                         ],
#     }
#     return source_enhancements

# def composite_defaults():
#     """Returns the default composite settings for each variable.

#     Returns
#     -------
#     dict
#         Dictionary in which keys are variable names and values are dictionaries
#         with `cmeta` settings.
#     """

#     cdefaults = {

#         'ndvi': {
#                     "composite_type": "mean",
#                     "temporal_interp": "linear",
#                     "spatial_interp": "nearest",
#                     "var_name": "ndvi",
#                     "var_unit": "-",
#                 },
#         'r0': {
#                     "composite_type": "mean",
#                     "temporal_interp": "linear",
#                     "spatial_interp": "nearest",
#                     "var_name": "r0",
#                     "var_unit": "-",
#                 },
#         'lulc': {
#                     "composite_type": 'mean',
#                     "temporal_interp": 'nearest',
#                     "spatial_interp": "nearest",
#                     "var_name": "lulc",
#                     "var_unit": "-",
#                 },
#         'z': {
#                     "composite_type": False,
#                     "temporal_interp": False,
#                     "spatial_interp": "bilinear",
#                     "var_name": "z",
#                     "var_unit": "m",
#                 },
#         'p': {
#                     "composite_type": "mean",
#                     "temporal_interp": "linear",
#                     "spatial_interp": "nearest",
#                     "var_name": "p_24",
#                     "var_unit": "mm/day",
#                 },
#         'se_root': {
#                     "composite_type": "max",
#                     "temporal_interp": False,
#                     "spatial_interp": "nearest",
#                     "var_name": "se_root",
#                     "var_unit": "-",
#                 },
#         'ra': {
#                     "composite_type": "mean",
#                     "temporal_interp": "linear",
#                     "spatial_interp": "bilinear",
#                     "var_name": "ra_24",
#                     "var_unit": "-",
#                 }
#         }

#     meteo = {
#             "composite_type": "mean",
#             "temporal_interp": "linear",
#             "spatial_interp": "bilinear",
#             }

#     meteo_names_units = {"t_air_max": "°C", "t_air_min":"°C", "t_air":"°C", 
#         "u2m":"m/s", "v2m":"m/s","qv":"kg/kg","p_air_0":"kPa"}

#     for name, unit in meteo_names_units.items():
#         cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **meteo}

#     static = {
#             "composite_type": False,
#             "temporal_interp": False,
#             "spatial_interp": "bilinear",    
#     }

#     statics_names_units = {'land_mask': "-",
#                         'lw_offset': "-",
#                         'lw_slope': "-",
#                         'r0_bare': "-",
#                         'r0_full': "-",
#                         'rn_offset': "-",
#                         'rn_slope': "-",
#                         'rs_min': "sm-1",
#                         't_amp_year': "C", 
#                         't_opt': "C",
#                         'vpd_slope': "mbar-1", 
#                         'z_obst_max': "m",
#                         'z_oro': "m",
#     }

#     for name, unit in statics_names_units.items():
#         cdefaults[name] = {**{"var_name": name, "var_unit": unit}, **static}

#     return cdefaults