def kelvin_to_celsius(ds, var, out_var = None):
    """Convert units of a variable from Kelvin to degrees Celsius.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the data.
    var : str
        Name of the variable inside ds to apply this function to.

    Returns
    -------
    xarray.Dataset
        Adjusted data.
    """

    new_data = ds[var] - 273.15
    new_data.attrs = {"unit": "°C"}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def celsius_to_kelvin(ds, var, out_var = None):
    """Convert units of a variable from Kelvin to degrees Celsius.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the data.
    var : str
        Name of the variable inside ds to apply this function to.

    Returns
    -------
    xarray.Dataset
        Adjusted data.
    """

    new_data = ds[var] + 273.15
    new_data.attrs = {"unit": "K"}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds