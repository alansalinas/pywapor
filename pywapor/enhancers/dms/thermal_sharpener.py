import os
import xarray as xr
from osgeo import gdal
import datetime
from pywapor.general.processing_functions import save_ds, remove_ds
from pywapor.enhancers.dms.pyDMS import DecisionTreeSharpener
import glob
import xml.etree.ElementTree as ET
import pandas as pd

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
        
    ds = xr.open_dataset(temporal_lr_input_data[0].split(":")[1])
    times = [x.strftime("%Y%m%d_%H%M%S") for x in pd.to_datetime(ds.time.values)]
    ds = ds.close()

    fps = list()
    for time_index, dt in enumerate(times):
        buildvrt_options = gdal.BuildVRTOptions(resolution = "highest", bandList = [time_index + 1])
        output_vrt = os.path.join(workdir, f"lowres_input_{dt}.vrt")
        if os.path.isfile(output_vrt):
            os.remove(output_vrt)
        ds = gdal.BuildVRT(output_vrt, temporal_lr_input_data, options = buildvrt_options)
        ds.FlushCache()
        ds = None
        fps.append(output_vrt)

    return fps

def preprocess(ds):
    dt = datetime.datetime.strptime(os.path.split(ds.encoding["source"])[-1], "highres_output_%Y%m%d_%H%M%S.nc")
    ds["Band1"] = ds["Band1"].expand_dims({"time": 1}).assign_coords({"time": [dt]})
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

    bt_nc = r"/Users/hmcoerver/Local/therm_sharpener/VIIRSL1/VNP02IMG.nc"
    folder = r"/Users/hmcoerver/Local/therm_sharpener"
    workdir = r"/Users/hmcoerver/Local/therm_sharpener/DMS"
    s2_ncs = glob.glob(os.path.join(folder, "*.nc"))
    z_nc = r"/Users/hmcoerver/Local/therm_sharpener/COPERNICUS/GLO90.nc"

    temporal_lr_input_data = f"NETCDF:{bt_nc}:bt"
    temporal_hr_input_data = [f"NETCDF:{fp}:{var}" for fp, var in [(fp, "_".join(os.path.split(fp)[-1].split("_")[0:-1])) for fp in s2_ncs]]
    static_hr_input_data = [f"NETCDF:{z_nc}:z", f"NETCDF:{z_nc}:slope", f"NETCDF:{z_nc}:aspect"]

    highResFiles = highres_inputs(workdir, temporal_hr_input_data, static_hr_input_data)
    lowResFiles = lowres_inputs(workdir, temporal_lr_input_data)

    fldr = r"/Users/hmcoerver/Local/therm_sharpener/DMS"

    var = "bt"

    highResFilenames = sorted(glob.glob(os.path.join(workdir, "highres*.vrt")))
    lowResFilenames = sorted(glob.glob(os.path.join(workdir, "lowres*.vrt")))

    out_fns = list()
    for highResFilename, lowResFilename in zip(highResFilenames, lowResFilenames):
        residualImage, correctedImage = thermal_sharpen(highResFilename, lowResFilename)
        
        fp, fn = os.path.split(highResFilename)
        fp_out = os.path.join(fp, fn.replace("input", "output").replace(".vrt", ".nc"))
        ds = gdal.Translate(fp_out, correctedImage)
        ds.FlushCache()
        ds = None

        out_fns.append(fp_out)

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