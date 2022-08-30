import pywapor
import datetime
from pywapor import et_look
import pywapor.general.levels as levels
from pywapor.general.logger import log, adjust_logger

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Local/20220325_20220415_test_data"

    # latlim = [28.9, 29.7]
    # lonlim = [30.2, 31.2]
    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]
    timelim = [datetime.date(2022, 4, 1), datetime.date(2022, 4, 11)]

    sources = "level_2_v3"
    sources = levels.pre_et_look_levels(sources)

    # sources = {k:v for k,v in sources.items() if k not in ["se_root", "r0"]}

    # sources = dict()
    # sources["ndvi"] = {'products': [{'source': 'SENTINEL2',
    #     'product_name': 'S2MSI2A',
    #     'enhancers': 'default',
    #     'is_example': True}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest',
    # 'composite_type': 'mean'}

    # sources["r0"] = {'products': [{'source': 'SENTINEL2',
    #     'product_name': 'S2MSI2A',
    #     'enhancers': 'default'}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest',
    # 'composite_type': 'min'}

    # sources["z"] = {'products': [{'source': 'COPERNICUS',
    #     'product_name': 'GLO30',
    #     'enhancers': 'default'}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest',
    # 'composite_type': 'max'}

    adjust_logger(True, folder, "INFO")

    ds = pywapor.pre_et_look.main(folder, latlim, lonlim, timelim, sources)
    assert ds.rio.crs.to_epsg() == 4326

    ds_out = pywapor.et_look.main(ds, et_look_version = "v3")
