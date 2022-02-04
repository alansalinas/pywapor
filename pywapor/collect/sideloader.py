import glob
import os
from pywapor.general.logger import log

def search_product_files(product_name, path):
    """Search for GeoTIFF files with a specific `product_name` in a folder
    defined by `path`.

    Parameters
    ----------
    product_name : str
        Should be inside a GeoTIFF file's name after one underscore.
    path : str
        Path to folder to search.

    Return
    ------
    list
        Paths to files with the product_name inside their filename.

    """
    fps = glob.glob(os.path.join(path, f"*_{product_name}_*.tif"))
    log.info(f"--> Collected {len(fps)} {product_name} file(s).")
    check = {os.path.split(fp)[-1]: check_filename(os.path.split(fp)[-1]) for fp in fps}
    log.debug(f"Filenames are incorrect ({check}).")
    return fps

def check_filename(fn):
    """Checks if files follow the correct naming conventions.         
    
    Filenames should look like this:
    `{variable_name}_{product_name}_{unit}_{period_length}_{date}.tif`
    where `variable_name`, `product_name` and `unit` are strings without 
    underscores. `period_length` can be either `-` (for instanteneous 
    data), an integer value or an integer value appended with "-daily". 
    `date` is either formatted as `%Y.%m.%d` in case a `period_length` is 
    specified OR `%Y.%m.%d.%H.%M` for instantenous data.

    Parameters
    ----------
    fn : str
        Filename to test.

    Returns
    -------
    bool
        True when the filename is correct.

    Raises
    ------
    ValueError
        The filename is incorrect.

    Examples
    --------
    >>> check_filenames("NDVI_MOD13_-_16-daily_2021.06.10")
    True
    >>> check_filenames("NDVI_MOD13_-_16_2021.06.10")
    True
    >>> check_filenames("ND_VI_MOD13_-_16_2021.06.10")
    False
    >>> check_filenames("NDVI_MOD13_-_16_2021.06.18.00.00")
    False
    >>> check_filenames("NDVI_MOD13_-_-_2021.06.18.00.00")
    True
    """
    name, ext = os.path.splitext(fn)
    parts = name.split("_")

    if ext != ".tif":
        raise ValueError("File has a wrong extension ({ext} != '.tif').")

    if len(parts) != 5:
        raise ValueError("Too many underscores in filename.")

    # valids = ["ALBEDO", "NDVI", "LST"]

    # if parts[0] not in valids:
    #     raise ValueError(f"Invalid variable name ({parts[0]}). Valid values are {valids}.")
    
    date_parts = parts[4].split(".")

    if len(date_parts) == 3:
        if "daily" in parts[3]:
            period_length = int(parts[3].split("-")[0])
            if period_length - float(parts[3].split("-")[0]) != 0:
                raise ValueError(f"Incorrect period_length ({period_length}).")
        else:
            if int(parts[3]) - float(parts[3]) != 0:
                raise ValueError(f"Incorrect period_length ({parts[3]}).")
    elif len(date_parts) == 5:
        if parts[3] != "-":
            raise ValueError(f"Incorrect period_length ({parts[3]}). Should be '-'.")
    else:
        raise ValueError(f"Invalid date ({parts[4]}).")

    return True
