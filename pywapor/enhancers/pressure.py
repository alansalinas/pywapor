def pa_to_kpa(ds, var):
    ds[var] = ds[var] / 1000
    ds[var].attrs = {}
    return ds