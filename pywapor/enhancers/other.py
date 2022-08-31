def drop_empty_times(ds, x, drop_vars = ["lst"], *args):
    for drop_var in drop_vars:
        ds = ds.isel(time = (ds[drop_var].count(dim=["y", "x"]) > 0).values)
    return ds