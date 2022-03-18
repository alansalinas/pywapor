import xarray as xr
import numpy as np
from pywapor.general.processing_functions import save_ds, domain_overlaps_domain

def opendap_to_xarray(store, fp, select, rename_keep_vars):
    
    # Open remote dataset.
    online_ds = xr.open_dataset(store)

    # Remove attributes.
    online_ds.attrs = {}

    # Check if selection is valid.
    checks = [domain_overlaps_domain(v, [online_ds[k].min(), online_ds[k].max()]) for k, v in select.items()]

    # Make the selection on the remote.
    if np.all(checks):
        online_ds = online_ds.sel({k: slice(*v) for k, v in select.items()})
    else:
        return None

    # Cleaning up.
    online_ds = online_ds[list(rename_keep_vars.keys())]
    online_ds = online_ds.rename_vars(rename_keep_vars)

    # Download the data.
    ds = save_ds(online_ds, fp)

    return ds