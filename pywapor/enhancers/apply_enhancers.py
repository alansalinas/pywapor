from functools import partial

def apply_enhancer(ds, variable, enhancer, source = None, log_it = True):
    ds = enhancer(ds, variable)

    if log_it:
        if isinstance(enhancer, partial):
            func_name = enhancer.func.__name__
            if "out_var" in enhancer.keywords:
                new_var = enhancer.keywords['out_var']
                label = f"--> Creating new variable `{new_var}`."
            else:
                label = f"--> Applying '{func_name}' to `{variable}`."
                if not isinstance(source, type(None)):
                    label = label[:-1] + f" from {source}."
        else:
            func_name = enhancer.__name__
            label = f"--> Applying '{func_name}' to `{variable}`."
            if not isinstance(source, type(None)):
                label = label[:-1] + f" from {source}."
    else:
        label = None

    return ds, label
