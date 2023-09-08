import os
import tempfile
import xarray as xr
from osgeo import gdal, osr
from pywapor.general.logger import log
from pywapor.general.processing_functions import save_ds, remove_ds
gdal.UseExceptions()

def calc_slope_or_aspect(ds, var, write_init = True):

    log.add()

    # Get filepath from xr.Dataset.
    fp = ds.encoding.get("source", None)

    if not var in ["slope", "aspect"]:
        raise ValueError("Please set `var` to either `'slope'` or `'aspect'`.")

    # Write data to file if necessary.
    if write_init or isinstance(fp, type(None)):
        fp = os.path.join(tempfile.gettempdir(), f"temp_{var}.nc")
        if os.path.isfile(fp):
            remove_ds(fp)
        while os.path.isfile(fp):
            fp = fp.replace(".nc", "_.nc")
        ds = save_ds(ds, fp, label = f"Generating {var} input file.")

    fp_out = f"/vsimem/temp_{var}.tif"

    ds_ = gdal.Open(fp)
    subdss_ = ds_.GetSubDatasets()

    if len(subdss_) > 0:
        fp = [x[0] for x in subdss_ if x[0].split(":")[-1] == "z"][0]
        ds_ = gdal.Open(fp)

    if ds_.RasterCount > 1:
        var_names = [ds_.GetRasterBand(i+1).GetMetadata().get("NETCDF_VARNAME", "unknown") for i in range(ds_.RasterCount)]
        band_select = var_names.find("z") if "z" in var_names else None
        if isinstance(band_select, type(None)):
            log.warning(f"--> Multiple bands found in `{fp}` and none is named `z`, using first band.")
            band_select = 1
    elif ds_.RasterCount == 1:
        band_select = 1
    else:
        raise ValueError(f"No bands found in `{fp}`.")

    band = ds_.GetRasterBand(band_select)
    wkt = ds_.GetProjection()
    src = osr.SpatialReference()
    src.ImportFromWkt(wkt)
    auth = src.GetAuthorityName(None)
    code = src.GetAuthorityCode(None)

    if f"{auth}:{code}" == "EPSG:4326":
        scale = band.GetScale()
        scale_ = {True: 1, False: scale}[isinstance(scale, type(None))] # NOTE in case there is no scale-factor, we set it to 1.
        scale_factor = int(111120 / scale_)
    else:
        scale_factor = 1
        log.warning("--> Dataset has different projection than `EPSG:4326`,\
                    slope scaling might be incorrect if CRS unit is not meters.")

    options = gdal.DEMProcessingOptions(
        slopeFormat = "degree",
        scale = scale_factor,
        computeEdges = True,
        zeroForFlat = True,
    )

    ds__ = gdal.DEMProcessing(fp_out, fp, var, options = options)
    ds__.FlushCache()
    ds__ = None

    ds_new = xr.open_dataset(fp_out).isel(band=0, drop=True).drop_vars(["x", "y"])
    ds_new = ds_new.rename({"band_data": var})
    ds_new.attrs = {}
    for var in ds_new.data_vars:
        ds_new[var].attrs = {}
    ds[var] = ds_new[var]

    log.sub()

    return ds
