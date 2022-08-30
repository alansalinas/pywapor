import numpy as np

def adjust_wind_height(ds, var, z = 10):
    """
    See FAO56 equation 47.
    """
    if "2m" in var:
        old_name = var.replace("2m", f"{z}m")
    else:
        old_name = var
    ds[var] = ds[old_name] * ((4.87)/(np.log(67.8*z - 5.42)))
    return ds

def windspeed(ds, var):
    if np.all([x in ds.data_vars for x in ["u10m", "v10m"]]):
        ds[var] = np.hypot(ds["u10m"], ds["v10m"])
    return ds

