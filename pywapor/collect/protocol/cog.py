import os
import tqdm
from pywapor.general.processing_functions import save_ds, open_ds, process_ds, remove_ds
from osgeo import gdal
import urllib
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.logger import log

def download(fp, product_name, coords, variables, post_processors, url_func, 
                gdal_config_options = {}, waitbar = True, ndv = -9999):

    folder, fn = os.path.split(fp)

    # Create folder.
    if not os.path.isdir(folder):
        os.makedirs(folder)

    for k, v in gdal_config_options.items():
        gdal.SetConfigOption(k, v)

    bands = [int(k.replace("Band", "")) for k in variables.keys() if "Band" in k]

    # Define bounding-box.
    bb = [coords["x"][1][0], coords["y"][1][1], 
            coords["x"][1][1], coords["y"][1][0]]

    options_dict = {
            "projWin": bb,
            "format": "netCDF",
            "creationOptions": ["COMPRESS=DEFLATE", "FORMAT=NC4"],
            "noData": ndv,
            "bandList": bands,
    }

    if waitbar:
        # Create waitbar.
        waitbar = tqdm.tqdm(position = 0, total = 100, bar_format='{l_bar}{bar}|', delay = 20)

        # Define callback function for waitbar progress.
        def _callback_func(info, *args):
            waitbar.update(info * 100 - waitbar.n)

        # Add callback to gdal.Translate options.
        options_dict.update({"callback": _callback_func})

    # Set gdal.Translate options.
    options = gdal.TranslateOptions(**options_dict)

    # Check if url is local or online path.
    url = url_func(product_name)
    is_not_local = urllib.parse.urlparse(url).scheme in ('http', 'https',)
    if is_not_local:
        url = f"/vsicurl/{url}"

    # Run gdal.Translate.
    ds = gdal.Translate(fp.replace(".nc", "_temp.nc"), url, options = options)

    # Reset the gdal.Dataset.
    ds.FlushCache()
    ds = None

    # Process the new netCDF.
    ds = open_ds(fp.replace(".nc", "_temp.nc"))

    ds = ds.rename_vars({k: f"Band{v}" for k,v in zip(list(ds.data_vars), bands)})

    ds = process_ds(ds, coords, variables)

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    # Save final output.
    out = save_ds(ds, fp, encoding = "initiate", label = f"Saving {fn}.")

    # Remove the temporary file.
    remove_ds(ds)

    return out
