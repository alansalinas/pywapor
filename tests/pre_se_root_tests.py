import pywapor
import datetime
from pywapor.general.logger import adjust_logger
import numpy as np

if __name__ == "__main__":

    ...
    ####
    ####

    sources = dict()
    sources["ndvi"] = {'products': [{'source': 'SENTINEL2',
        'product_name': 'S2MSI2A',
        'enhancers': 'default',
        'is_example': True}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest'}

    sources["lst"] = {'products': [{'source': 'SENTINEL3',
        'product_name': 'SL_2_LST___',
        'enhancers': 'default'}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest'}

    sources["bt"] = {'products': [{'source': 'VIIRSL1',
        'product_name': 'VNP02IMG',
        'enhancers': 'default'}],
    'temporal_interp': 'linear',
    'spatial_interp': 'nearest'}

    folder = r"/Users/hmcoerver/Local/test1"
    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]
    timelim = [datetime.date(2022, 4, 1), datetime.date(2022, 4, 11)]
    adjust_logger(True, folder, "INFO")
    ds = pywapor.pre_se_root.main(folder, latlim, lonlim, timelim, sources)
    assert ds.rio.crs.to_epsg() == 4326
    assert "bt" not in ds.data_vars
    assert "lst" in ds.data_vars
    assert "ndvi" in ds.data_vars
    assert ds.ndvi.min().values >= -1.
    assert ds.ndvi.max().values <= 1.
    assert 240 < ds.lst.mean().values < 330

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test2"
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2022, 4, 1), datetime.date(2022, 4, 3)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_2_v3"
    input_data = pywapor.pre_se_root.main(folder, latlim, lonlim, timelim, 
                                            sources, bin_length = 3)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "p_air_i", "p_air_0_i", "r0_bare", "r0_full", "t_air_i", "t_dew_i", "u_i", "wv_i", "lst"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert 90 < input_data.p_air_i.mean().values < 110
    assert 90 < input_data.p_air_0_i.mean().values < 110
    assert 0 < input_data.r0_bare.mean().values < 1
    assert 0 < input_data.r0_full.mean().values < 1
    assert -40 < input_data.t_air_i.mean().values < 50
    assert -40 < input_data.t_dew_i.mean().values < 50
    assert 0 < input_data.u_i.mean().values < 150
    assert 0 < input_data.wv_i.mean().values < 100
    assert 240 < input_data.lst.mean().values < 320

    ds = pywapor.se_root.main(input_data, se_root_version = "v3")
    assert ds.rio.crs.to_epsg() == 4326
    assert "se_root" in ds.data_vars
    assert ds.se_root.min().values >= 0.
    assert ds.se_root.max().values <= 1.

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test3"
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2019, 4, 1), datetime.date(2019, 4, 11)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_1"
    input_data = pywapor.pre_se_root.main(folder, latlim, lonlim, timelim, sources)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "lst"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert 90 < input_data.p_air_i.mean().values < 110
    assert 90 < input_data.p_air_0_i.mean().values < 110
    assert 0 < input_data.r0_bare.mean().values < 1
    assert 0 < input_data.r0_full.mean().values < 1
    assert -40 < input_data.t_air_i.mean().values < 50
    assert 0 < input_data.wv_i.mean().values < 100
    assert 240 < input_data.lst.mean().values < 320
    
    ds = pywapor.se_root.main(input_data, se_root_version = "v2")
    assert ds.rio.crs.to_epsg() == 4326
    assert "se_root" in ds.data_vars
    assert ds.se_root.min().values >= 0.
    assert ds.se_root.max().values <= 1.

    ####
    ####

    folder = r"/Users/hmcoerver/Local/test4"
    adjust_logger(True, folder, "INFO")
    timelim = [datetime.date(2019, 4, 1), datetime.date(2019, 4, 3)]
    latlim = [29.4, 29.6]
    lonlim = [30.7, 30.9]
    sources = "level_2"
    input_data = pywapor.pre_se_root.main(folder, latlim, lonlim, timelim, 
                                            sources, bin_length = 3)
    assert input_data.rio.crs.to_epsg() == 4326
    assert np.all([x in input_data.data_vars for x in ["ndvi", "lst"]])
    assert input_data.ndvi.min().values >= -1.
    assert input_data.ndvi.max().values <= 1.
    assert 90 < input_data.p_air_i.mean().values < 110
    assert 90 < input_data.p_air_0_i.mean().values < 110
    assert 0 < input_data.r0_bare.mean().values < 1
    assert 0 < input_data.r0_full.mean().values < 1
    assert -40 < input_data.t_air_i.mean().values < 50
    assert 0 < input_data.wv_i.mean().values < 100
    assert 240 < input_data.lst.mean().values < 320

    ds = pywapor.se_root.main(input_data, se_root_version = "v2")
    assert ds.rio.crs.to_epsg() == 4326
    assert "se_root" in ds.data_vars
    assert ds.se_root.min().values >= 0.
    assert ds.se_root.max().values <= 1.
