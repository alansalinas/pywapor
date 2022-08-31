import pywapor
import datetime
from pywapor import et_look
import pywapor.general.levels as levels
from pywapor.general.logger import log, adjust_logger
import numpy as np

if __name__ == "__main__":

    sources = dict()
    sources["ndvi"] = {'products': [{'source': 'SENTINEL2',
        'product_name': 'S2MSI2A',
        'enhancers': 'default',
        'is_example': True}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest',
    'composite_type': 'mean'}

    sources["r0"] = {'products': [{'source': 'SENTINEL2',
        'product_name': 'S2MSI2A',
        'enhancers': 'default'}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest',
    'composite_type': 'min'}

    sources["z"] = {'products': [{'source': 'COPERNICUS',
        'product_name': 'GLO30',
        'enhancers': 'default'}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest',
    'composite_type': 'max'}

    folder = r"/Users/hmcoerver/Local/test5" # ✅
    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]
    timelim = [datetime.date(2022, 4, 1), datetime.date(2022, 4, 11)]
    adjust_logger(True, folder, "INFO")
    ds = pywapor.pre_et_look.main(folder, latlim, lonlim, timelim, sources)
    assert ds.rio.crs.to_epsg() == 4326
    assert "ndvi" in ds.data_vars
    assert "r0" in ds.data_vars
    assert "ndvi" in ds.data_vars
    assert ds.ndvi.min().values >= -1.
    assert ds.ndvi.max().values <= 1.

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test6" # ✅
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2022, 4, 5), datetime.date(2022, 4, 9)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_2_v3"
    input_data = pywapor.pre_et_look.main(folder, latlim, lonlim, timelim, 
                                            sources, bin_length = 4)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "r0", "se_root"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert input_data.r0.min().values >= 0.
    assert input_data.r0.max().values <= 1.
    assert input_data.se_root.min().values >= 0.
    assert input_data.se_root.max().values <= 1.

    ds = pywapor.et_look.main(input_data, et_look_version = "v3")
    assert ds.rio.crs.to_epsg() == 4326
    assert "et_24_mm" in ds.data_vars
    assert ds.et_24_mm.min().values >= 0.
    assert ds.et_24_mm.max().values <= 16.

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test7" # ✅
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2021, 7, 1), datetime.date(2021, 7, 11)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_1"
    input_data = pywapor.pre_et_look.main(folder, latlim, lonlim, timelim, 
                                            sources)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "r0", "se_root"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert input_data.r0.min().values >= 0.
    assert input_data.r0.max().values <= 1.
    assert input_data.se_root.min().values >= 0.
    assert input_data.se_root.max().values <= 1.

    ds = pywapor.et_look.main(input_data, et_look_version = "v2")
    assert ds.rio.crs.to_epsg() == 4326
    assert "et_24_mm" in ds.data_vars
    assert ds.et_24_mm.min().values >= 0.
    assert ds.et_24_mm.max().values <= 16.

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test8" # ✅
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2019, 4, 1), datetime.date(2019, 4, 11)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_2"
    input_data = pywapor.pre_et_look.main(folder, latlim, lonlim, timelim, 
                                            sources)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "r0", "se_root"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert input_data.r0.min().values >= 0.
    assert input_data.r0.max().values <= 1.
    assert input_data.se_root.min().values >= 0.
    assert input_data.se_root.max().values <= 1.

    ds = pywapor.et_look.main(input_data, et_look_version = "v2")
    assert ds.rio.crs.to_epsg() == 4326
    assert "et_24_mm" in ds.data_vars
    assert ds.et_24_mm.min().values >= 0.
    assert ds.et_24_mm.max().values <= 16.
