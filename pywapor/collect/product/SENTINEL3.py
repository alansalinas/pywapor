def default_vars(product_name, req_vars):

    variables = {

        "SL_2_LST___": {
                    "LST_in.nc": [(), "lst"],
                    "geodetic_in.nc": [(), "coords"],
                    "flags_in.nc": [(), "qa"],
                    "time_in.nc": [(), "time"]
        }
        
    }

    req_dl_vars = {

        "SL_2_LST___": {
            "lst": ["LST_in.nc", "geodetic_in.nc", "flags_in.nc", "time_in.nc"],
        },
    }

    out = {val:variables[product_name][val] for sublist in map(req_dl_vars[product_name].get, req_vars) for val in sublist}

    return out