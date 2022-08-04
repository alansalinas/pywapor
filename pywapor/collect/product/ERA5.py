import os
from pywapor.collect.protocol import cds
from pywapor.general.processing_functions import open_ds
from pywapor.general.logger import log
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.enhancers.pressure import pa_to_kpa
from pywapor.enhancers.temperature import kelvin_to_celsius

def default_vars(product_name, req_vars):

    variables = {
        "sis-agrometeorological-indicators": {
                    "2m_temperature": [{"statistic": "24_hour_mean", "format": "zip"}, "t_air"],
                    "2m_dewpoint_temperature": [{"statistic": "24_hour_mean", "format": "zip"}, "t_dew"],
                    "2m_relative_humidity": [{"time": [f"{x:02d}_00" for x in [6,9,12,15,18]], "format": "zip"}, "rh"],
                    "10m_wind_speed": [{"statistic": "24_hour_mean", "format": "zip"}, "u"],
                    "vapour_pressure": [{"statistic": "24_hour_mean", "format": "zip"}, "vp"],
                    "solar_radiation_flux": [{"statistic": "24_hour_mean", "format": "zip"}, "ra"],
                        },

        "reanalysis-era5-single-levels": {
            '10m_u_component_of_wind':  [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "u_10m"],
            '10m_v_component_of_wind':  [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "v_10m"],
            '2m_dewpoint_temperature':  [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "t_dew"],
            'mean_sea_level_pressure':  [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "p_air_0"],
            'surface_pressure':         [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "p_air"],
            '2m_temperature':           [{"time": [f"{x:02d}:00" for x in range(24)], "format": "netcdf", "product_type": "reanalysis"}, "t_air"],
        }
    }

    req_dl_vars = {
        "sis-agrometeorological-indicators": {
            "t_air": ["2m_temperature"],
            "t_dew": ["2m_dewpoint_temperature"],
            "rh": ["2m_relative_humidity"],
            "u": ["10m_wind_speed"],
            "vp": ["vapour_pressure"],
            "ra": ["solar_radiation_flux"],
        },
        "reanalysis-era5-single-levels": {
            "u_10m": ['10m_u_component_of_wind'],
            "v_10m": ['10m_v_component_of_wind'],
            "t_dew": ['2m_dewpoint_temperature'],
            "p_air_0": ['mean_sea_level_pressure'],
            "p_air": ['surface_pressure'],
            "t_air": ['2m_temperature'],
        }
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out

def default_post_processors(product_name, req_vars):

    # TODO check if all units are correct (!!!)
    post_processors = {
        "sis-agrometeorological-indicators": {
            "t_air": [kelvin_to_celsius],
            "t_dew": [],
            "rh": [],
            "u": [], # TODO convert 10m to 2m (!!!)
            "vp": [],
            "ra": [],
        },
        "reanalysis-era5-single-levels": {
            "u_10m": [], # TODO convert 10m to 2m (!!!)
            "v_10m": [], # TODO convert 10m to 2m (!!!)
            "t_dew": [],
            "p_air_0": [pa_to_kpa],
            "p_air": [pa_to_kpa],
            "t_air": [kelvin_to_celsius],
        }
    }

    out = {k:v for k,v in post_processors[product_name].items() if k in req_vars}

    return out

def download(folder, latlim, lonlim, timelim, product_name, req_vars, 
                variables = None, post_processors = None):

    product_folder = os.path.join(folder, "ERA5")

    if not os.path.exists(product_folder):
        os.makedirs(product_folder)

    fn_final = os.path.join(product_folder, f"{product_name}.nc")
    if os.path.isfile(fn_final):
        return open_ds(fn_final, "all")

    spatial_buffer = True
    if spatial_buffer:
        latlim = [latlim[0] - 0.1, latlim[1] + 0.1]
        lonlim = [lonlim[0] - 0.1, lonlim[1] + 0.1]

    if isinstance(variables, type(None)):
        variables = default_vars(product_name, req_vars)

    if isinstance(post_processors, type(None)):
        post_processors = default_post_processors(product_name, req_vars)
    else:
        default_processors = default_post_processors(product_name, req_vars)
        post_processors = {k: {True: default_processors[k], False: v}[v == "default"] for k,v in post_processors.items()}

    ds = cds.download(product_folder, product_name, latlim, lonlim, timelim, variables)

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    return ds

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/On My Mac/era_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = ["2021-06-26", "2021-07-05"]

    product_name = "sis-agrometeorological-indicators"
    # product_name = "reanalysis-era5-single-levels"
    req_vars = ["t_air", "t_dew", "rh", "u", "vp", "ra"]
    # req_vars = ["u_10m", "v_10m", "t_dew", "p_air_0", "p_air", "t_air"]

    variables = None
    post_processors = None

    ds = download(folder, latlim, lonlim, timelim, product_name = product_name, 
                req_vars = req_vars, variables = variables, post_processors = post_processors)
