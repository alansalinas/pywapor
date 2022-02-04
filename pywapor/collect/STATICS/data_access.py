from osgeo import gdal
import tqdm
import os
import urllib
from pywapor.general.logger import log
from pywapor.general.pre_defaults import composite_defaults

def get_cog(url, out_file, bb, srs = "EPSG:4326"):

    # Check if out_file folder exists.
    folder = os.path.split(out_file)[0]
    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create waitbar
    waitbar = tqdm.tqdm(desc= f"Tile: 0 / 1",
                                position = 0,
                                total = 100,
                                bar_format='{l_bar}{bar}|')

    # Advance tile counter.
    waitbar_i = int(waitbar.desc.split(" ")[1])
    waitbar_desc = str(waitbar.desc)
    waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /")) 

    # Define callback function for waitbar progress.
    def callback_func(info, *args):
        waitbar.update(info * 100 - waitbar.n)

    # Set gdal.Translate options.
    options = gdal.TranslateOptions(
        projWin = bb,
        projWinSRS = srs,
        outputSRS = srs,
        callback = callback_func,
    )

    # Check if url is local or online path.
    is_not_local = urllib.parse.urlparse(url).scheme in ('http', 'https',)
    if is_not_local:
        url = f"/vsicurl/{url}"

    # Run gdal.Translate.
    ds = gdal.Translate(out_file, url, options = options)

    # Finish waitbar.
    waitbar.refresh()

    # Reset the gdal.Dataset.
    ds.FlushCache()
    ds = None

    return out_file

def collect(Dir, latlim, lonlim, vars = ['land_mask', 'lw_offset',
            'lw_slope', 'r0_bare', 'r0_full', 'rn_offset', 'rn_slope',
            'rs_min', 't_amp_year', 't_opt', 'vpd_slope', 'z_obst_max',
            'z_oro'], base_url = r"https://storage.googleapis.com/fao-cog-data", **kwargs):

    bb = [lonlim[0], latlim[1], lonlim[1], latlim[0]]

    files = list()

    log("--> Downloading STATICS.")

    cmetas = composite_defaults()

    for var in vars:

        var_name = var.replace("_", "-")

        if var in cmetas.keys():
            unit = cmetas[var]["var_unit"]
            fn = f"{var_name}_STATICS_{unit}_-_-.tif"
        else:
            fn = f"{var_name}_STATICS_unknown_-_-.tif"

        out_file = os.path.join(Dir, "STATICS", fn)
        url = os.path.join(base_url, f"L1_{var}.cog.tif")

        if not os.path.isfile(out_file):
            out_file = get_cog(url, out_file, bb)

        files.append(out_file)

    return files

if __name__ == "__main__":

    # url = r"http://landsat-pds.s3.amazonaws.com/c1/L8/073/087/LC08_L1TP_073087_20190818_20190902_01_T1/LC08_L1TP_073087_20190818_20190902_01_T1_B2.TIF"
    # url = r"/Volumes/Data/L1_STATICS_cogs/L1_land_mask.cog.tif"
    # bb = [175.2, -38.8, 175.8, -39.5]

    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]

    raw_folder = r"/Users/hmcoerver/On My Mac/cog_test"

    out = collect(raw_folder, latlim, lonlim)
