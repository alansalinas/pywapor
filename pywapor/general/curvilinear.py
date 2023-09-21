import numpy as np
import os
from osgeo import gdal
from collections import OrderedDict
from pywapor.general.logger import log

def curvi_to_recto(lats, lons, data, out_fn, warp_kwargs = {}):

    if isinstance(lats, str):
        lats = gdal.Open(lats)
    if isinstance(lons, str):
        lons = gdal.Open(lons)

    for name, x in data.items():
        if isinstance(x, str):
            data[name] = gdal.Open(x)

    if lats.RasterCount != 1:
        raise ValueError
    if lons.RasterCount != 1:
        raise ValueError
    
    xsize = lons.RasterXSize
    ysize = lons.RasterYSize
    if xsize != lats.RasterXSize:
        raise ValueError
    if ysize != lons.RasterYSize:
        raise ValueError

    data_filtered = OrderedDict()
    for name, x in data.items():
        if x.RasterCount > 0:
            data_filtered[name] = x
        else:
            log.warning(f"No Raster found in {name} ({x.GetFileList()})")

    #####
    # Create VRT file linking all the required data.
    #####
    vrt_options = gdal.BuildVRTOptions(separate = True)
    data_bands = gdal.BuildVRT("/vsimem/temp.vrt", list(data_filtered.values()), options = vrt_options) 

    items = {"LINE_OFFSET": "0", 
            "LINE_STEP": "1", 
            "PIXEL_OFFSET": "0", 
            "PIXEL_STEP": "1", 
            "X_BAND": "1",
            "X_DATASET": lons.GetDescription(),
            "Y_BAND": "1",
            "Y_DATASET": lats.GetDescription()}

    _ = data_bands.SetMetadata(items, "GEOLOCATION")
    data_bands.FlushCache()

    #####
    # Warp the data.
    #####
    coptions = [
        "COMPRESS=DEFLATE"
        ]

    options = gdal.WarpOptions(
        srcSRS = "epsg:4326",
        dstSRS = "epsg:4326",
        dstNodata = -9999,
        geoloc = True,
        creationOptions = coptions,
        **warp_kwargs,
    )

    try:
        # NOTE turn of logging because of irrelevant ERROR message. "NUL" for windows, other for mac/linux.
        # also see https://gis.stackexchange.com/questions/358304/suppress-warnings-in-python-ogr
        gdal.SetConfigOption("CPL_LOG", {"nt": "NUL"}.get(os.name, "/dev/null"))
        out = gdal.Warp(out_fn, data_bands, options = options)
    except Exception as e:
        raise e
    finally:
        gdal.SetConfigOption("CPL_LOG", None)

    _ = [out.GetRasterBand(i + 1).SetMetadata({"variable": list(data_filtered.keys())[i]}) for i in range(out.RasterCount)]
    out.FlushCache()
    out = None

    return out_fn

def create_grid(latlim, lonlim, dx_dy = (0.0033, 0.0033)):
    dx, dy = dx_dy
    nx = np.ceil((lonlim[1] - lonlim[0]) / dx)
    ny = np.ceil((latlim[1] - latlim[0]) / dy)
    bb = [lonlim[0], latlim[0], lonlim[0] + nx * dx, latlim[0] + ny * dy]
    return bb, nx, ny
