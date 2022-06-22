import os
import tqdm
from pywapor.general.processing_functions import save_ds, open_ds, process_ds
from osgeo import gdal
import urllib
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.general.logger import log

def download(folder, product_name, coords, variables, post_processors, url_func):

    # Create folder.
    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create waitbar.
    waitbar = tqdm.tqdm(position = 0, total = 100, bar_format='{l_bar}{bar}|', delay = 20)

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
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)
    
    # Save final output.
    ds = save_ds(ds, fp, decode_coords = "all")

    # Remove the temporary file.
    os.remove(fp.replace(".nc", "_temp.nc"))

    return ds

# if __name__ == "__main__":

#     folder = r"/Users/hmcoerver/Downloads/pywapor_test/GLOBCOVER"
#     product_name = r"GLOBCOVER"

#     latlim = [28.9, 29.7]
#     lonlim = [30.2, 31.2]

#     coords = {"x": ("lon", lonlim), "y": ("lat", latlim)}

#     variables = {
#                 "Band1": [("lat", "lon"), "lulc"],
#                 "crs": [(), "spatial_ref"],
#                     }


#     post_processors = []

#     def url_func(product_name):
#         return r"http://due.esrin.esa.int/files/GLOBCOVER_L4_200901_200912_V2.3.color.tif"

    # download(folder, product_name, coords, variables, post_processors, url_func)