from functools import partial

def apply_enhancer(ds, variable, enhancer, source = None, log_it = True):
    """_summary_

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to enhance.
    variable : str
        Variable name to enhance.
    enhancer : function
        Function to be applied to `variable` inside `ds`. Should take `ds` as
        first argument and `variable` as second.
    source : str, optional
        Specifies the source of the data inside `ds` for logging purposes, by default None.
    log_it : bool, optional
        Create a label to log or not, by default True.

    Returns
    -------
    xr.Dataset
        The enhanced dataset.
    str
        The label to log when calculating the dataset.
    """
    ds = enhancer(ds, variable)

    if log_it:
        if isinstance(enhancer, partial):
            func_name = enhancer.func.__name__
            if "out_var" in enhancer.keywords:
                new_var = enhancer.keywords['out_var']
                label = f"--> Creating new variable `{new_var}`."
            else:
                label = f"--> Applying '{func_name}'."
                if not isinstance(source, type(None)):
                    label = label[:-1] + f" from {source}."
        else:
            func_name = enhancer.__name__
            label = f"--> Applying '{func_name}'."
            if not isinstance(source, type(None)):
                label = label[:-1] + f" from {source}."
    else:
        label = None

    return ds, label
