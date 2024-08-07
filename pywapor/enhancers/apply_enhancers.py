from functools import partial
import functools
from pywapor.general.logger import log
from pywapor.general.processing_functions import func_from_string

def apply_enhancers(post_processors, ds):
    for var, funcs in post_processors.items():
        for func in funcs:
            if isinstance(func, str):
                func = func_from_string(func)
            elif isinstance(func, dict):
                args = func.get("args", ())
                kwargs = func.get("keywords", {})
                func = func_from_string(func["func"])
                if kwargs or args:
                    func = functools.partial(func, *args, **kwargs)
            ds, label = apply_enhancer(ds, var, func)
    return ds

def apply_enhancer(ds, variable, enhancer):
    """Apple a function to a (variable in a) dataset. 

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to enhance.
    variable : str
        Variable name to enhance.
    enhancer : function
        Function to be applied to `variable` inside `ds`. Should take `ds` as
        first argument and `variable` as second.

    Returns
    -------
    tuple
        The enhanced dataset and the label to log when calculating the dataset.
    """
    if isinstance(enhancer, partial):
        func_name = enhancer.func.__name__
    else:
        func_name = enhancer.__name__
    
    if isinstance(variable, type(None)):
        label = f"--> Applying '{func_name}'."
    else:
        label = f"--> Applying '{func_name}' to `{variable}`."

    log.info(label)

    ds = enhancer(ds, variable)

    return ds, label
