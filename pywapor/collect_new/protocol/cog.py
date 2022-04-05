import os
import tqdm
from pywapor.general.processing_functions import save_ds, open_ds, process_ds
from osgeo import gdal
import urllib

def download(folder, product_name, coords, variables, post_processors, url_func):

    # Create waitbar
    waitbar = tqdm.tqdm(position = 0, total = 100, bar_format='{l_bar}{bar}|')

    # Define callback function for waitbar progress.
    def _callback_func(info, *args):
        waitbar.update(info * 100 - waitbar.n)

    # Define filepath.
    fp = os.path.join(folder, f"{product_name}.nc")

    # Define bounding-box.
    bb = [coords["x"][1][0], coords["y"][1][1], 
            coords["x"][1][1], coords["y"][1][0]]

    # Set gdal.Translate options.
    options = gdal.TranslateOptions(
        projWin = bb,
        format = "netCDF",
        creationOptions = ["COMPRESS=DEFLATE", "FORMAT=NC4C"],
        callback = _callback_func,
    )

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
    ds = open_ds(fp.replace(".nc", "_temp.nc"), decode_coords = "all")
    ds = process_ds(ds, coords, variables)

    # Apply product specific functions.
    for func in post_processors:
        ds = func(ds)
    
    # Save final output.
    ds = save_ds(ds, fp, decode_coords = "all")

    # Remove the temporary file.
    os.remove(fp.replace(".nc", "_temp.nc"))

    return ds
