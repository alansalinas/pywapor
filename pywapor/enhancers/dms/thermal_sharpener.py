import os
import xarray as xr
from osgeo import gdal
import datetime
from pywapor.general.processing_functions import save_ds, remove_ds, open_ds
from pywapor.enhancers.dms.pyDMS import DecisionTreeSharpener
import xml.etree.ElementTree as ET
import pandas as pd
from pywapor.general.logger import log, adjust_logger

def gdal_stop_warnings_and_raise_errors():
    """https://gis.stackexchange.com/questions/43404/how-to-detect-a-gdal-ogr-warning/68042
    """
    class GdalErrorHandler(object):
        def __init__(self):
            self.err_level=gdal.CE_None
            self.err_no=0
            self.err_msg=''
        def handler(self, err_level, err_no, err_msg):
            self.err_level=err_level
            self.err_no=err_no
            self.err_msg=err_msg
    err=GdalErrorHandler()
    handler=err.handler
    gdal.PushErrorHandler(handler)
    gdal.UseExceptions()

gdal_stop_warnings_and_raise_errors()

def highres_inputs(workdir, temporal_hr_input_data, static_hr_input_data):

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    buildvrt_options = gdal.BuildVRTOptions(resolution = "highest", separate = True)
    output_vrt = os.path.join(workdir, "highres_input_template.vrt")
    if os.path.isfile(output_vrt):
        os.remove(output_vrt)
    ds = gdal.BuildVRT(output_vrt, temporal_hr_input_data + static_hr_input_data, options = buildvrt_options)
    ds.FlushCache()
    ds = None

    ds = xr.open_dataset(temporal_hr_input_data[0].split(":")[1])
    times = [x.strftime("%Y%m%d_%H%M%S") for x in pd.to_datetime(ds.time.values)]
    ds = ds.close()

    with open(output_vrt, 'r') as f:
        data = f.read()

    fps = list()
    for time_index, dt in enumerate(times):
        tree = ET.fromstring(data)
        for x in tree.findall("VRTRasterBand"):
            source_node = x.findall("ComplexSource")[0]
            fn_node = source_node.findall("SourceFilename")[0]
            if fn_node.text in temporal_hr_input_data:
                band_node = source_node.findall("SourceBand")[0]
                band_node.text = str(time_index + 1)

        fp_out = output_vrt.replace("_template.vrt", f"_{dt}.vrt")
        with open(fp_out, 'wb') as output:
            output.write(ET.tostring(tree))
        fps.append(fp_out)

    os.remove(output_vrt)

    return fps

def lowres_inputs(workdir, temporal_lr_input_data):

    if not os.path.exists(workdir):
        os.makedirs(workdir)
        
    ds = xr.open_dataset(temporal_lr_input_data.split(":")[1])
    times = [x.strftime("%Y%m%d_%H%M%S") for x in pd.to_datetime(ds.time.values)]
    ds = ds.close()

    fps = list()
    for time_index, dt in enumerate(times):
        buildvrt_options = gdal.BuildVRTOptions(resolution = "highest", bandList = [time_index + 1])
        output_vrt = os.path.join(workdir, f"lowres_input_{dt}.vrt")
        if os.path.isfile(output_vrt):
            os.remove(output_vrt)
        ds = gdal.BuildVRT(output_vrt, [temporal_lr_input_data], options = buildvrt_options)
        ds.FlushCache()
        ds = None
        fps.append(output_vrt)

    return fps

def preprocess(ds):
    dt = datetime.datetime.strptime(os.path.split(ds.encoding["source"])[-1], "highres_output_%Y%m%d_%H%M%S.nc")
    ds["Band1"] = ds["Band1"].expand_dims({"time": 1}).assign_coords({"time": [dt]})
    return ds

def sharpen(ds, var, dss, folder, vars_for_sharpening = ['nmdi', 'bsi', 'mndwi', 
                'vari_red_edge', 'psri', 'nir', 'green', 'z', 'aspect', 'slope']):

    missing_vars = [x for x in vars_for_sharpening if x not in dss.keys()]
    if len(missing_vars) > 0:
        log.info(f"--> Sharpening `{var}` without `{'` and `'.join(missing_vars)}`.")

    workdir = os.path.join(folder, "DMS")

    temporal_lr_input_data = f"NETCDF:{ds.encoding['source']}:{var}"
    temporal_hr_input_data = [f"NETCDF:{ds.encoding['source']}:{x}" for x, ds in dss.items() if "time" in ds.dims and x in vars_for_sharpening]
    static_hr_input_data = [f"NETCDF:{ds.encoding['source']}:{x}" for x, ds in dss.items() if "time" not in ds.dims and x in vars_for_sharpening]

    highResFiles = sorted(highres_inputs(workdir, temporal_hr_input_data, static_hr_input_data))
    lowResFiles = sorted(lowres_inputs(workdir, temporal_lr_input_data))

    out_fns = list()
    log.info(f"--> Sharpening {len(highResFiles)} `{var}` images.").add()
    for i, (highResFilename, lowResFilename) in enumerate(zip(highResFiles, lowResFiles)):
        fp, fn = os.path.split(highResFilename)

        log.info(f"--> ({i+1}/{len(highResFiles)}) Sharpening `{os.path.split(lowResFilename)[-1]}` with `{fn}`.").add()

        _, correctedImage = thermal_sharpen(highResFilename, lowResFilename)

        log.sub()

        fp_out = os.path.join(fp, fn.replace("input", "output").replace(".vrt", ".nc"))
        ds = gdal.Translate(fp_out, correctedImage)
        ds.FlushCache()
        ds = None

        out_fns.append(fp_out)

    log.sub()

    ds = xr.open_mfdataset(out_fns, preprocess = preprocess, decode_coords = "all")
    ds = ds.rename_dims({"lat": "y", "lon": "x"})
    ds = ds.transpose("time", "y", "x")
    ds = ds.rename_vars({"crs": "spatial_ref", "Band1": var, "lat": "y", "lon": "x"})
    ds = ds.rio.write_grid_mapping("spatial_ref")
    ds = ds.sortby("y", ascending = False)
    ds = ds.sortby("x")
    ds.attrs = {}

    fp = os.path.join(workdir, f"{var}_i.nc")

    ds = save_ds(ds, fp, encoding = "initiate", label = f"Merging sharpened `{var}` files.")

    for fp in out_fns:
        remove_ds(fp)

    return ds

def thermal_sharpen(highres_fn, lowres_fn):

    commonOpts = {
                "highResFiles":                     [highres_fn],
                "lowResFiles":                      [lowres_fn],
                "lowResQualityFiles":               [],
                "lowResGoodQualityFlags":           [],
                "cvHomogeneityThreshold":           0.20,
                "movingWindowSize":                 0,
                "disaggregatingTemperature":        True
            }

    dtOpts = {
            "perLeafLinearRegression":              True,
            "linearRegressionExtrapolationRatio":   0.25
            }

    opts = commonOpts.copy()
    opts.update(dtOpts)
    disaggregator = DecisionTreeSharpener(**opts)

    disaggregator.trainSharpener()
    downscaledFile = disaggregator.applySharpener(highres_fn, lowres_fn)
    residualImage, correctedImage = disaggregator.residualAnalysis(downscaledFile, lowres_fn)

    return residualImage, correctedImage

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Local/therm_sharpener"

    adjust_logger(True, folder, "INFO")

    vars_for_sharpening = ['nmdi', 'bsi', 'mndwi', 'vari_red_edge', 'psri', 'nir', 'green', 'z', 'aspect', 'slope']

    var = "bt"
    ds = open_ds(r"/Users/hmcoerver/Local/therm_sharpener/VIIRSL1/VNP02IMG.nc")

    dss = {
        'nmdi':             open_ds('/Users/hmcoerver/Local/therm_sharpener/nmdi_i.nc'),
        'bsi':              open_ds('/Users/hmcoerver/Local/therm_sharpener/bsi_i.nc'),
        'mndwi':            open_ds('/Users/hmcoerver/Local/therm_sharpener/mndwi_i.nc'),
        'vari_red_edge':    open_ds('/Users/hmcoerver/Local/therm_sharpener/vari_red_edge_i.nc'),
        'psri':             open_ds('/Users/hmcoerver/Local/therm_sharpener/psri_i.nc'),
        'nir':              open_ds('/Users/hmcoerver/Local/therm_sharpener/nir_i.nc'),
        'green':            open_ds('/Users/hmcoerver/Local/therm_sharpener/green_i.nc'),
        'z':                open_ds('/Users/hmcoerver/Local/therm_sharpener/COPERNICUS/GLO90.nc'),
        'aspect':           open_ds('/Users/hmcoerver/Local/therm_sharpener/COPERNICUS/GLO90.nc'),
        'slope':            open_ds('/Users/hmcoerver/Local/therm_sharpener/COPERNICUS/GLO90.nc'),
        }

    out = sharpen(ds, var, dss, folder)