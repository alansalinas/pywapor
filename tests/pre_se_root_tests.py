import pywapor
import datetime
import pywapor.general.levels as levels
from pywapor.general.logger import log, adjust_logger

if __name__ == "__main__":

    folder = r"/Users/hmcoerver/Local/20220325_20220415_test_data"

    latlim = [29.4, 29.7]
    lonlim = [30.7, 31.0]
    timelim = [datetime.date(2022, 4, 1), datetime.date(2022, 4, 11)]

    sources = "level_2_v3"
    sources = levels.pre_se_root_levels(sources)

    # sources = dict()
    # sources["ndvi"] = {'products': [{'source': 'SENTINEL2',
    #     'product_name': 'S2MSI2A',
    #     'enhancers': 'default',
    #     'is_example': True}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest'}

    # sources["lst"] = {'products': [{'source': 'SENTINEL3',
    #     'product_name': 'SL_2_LST___',
    #     'enhancers': 'default'}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest'}

    # sources["bt"] = {'products': [{'source': 'VIIRSL1',
    #     'product_name': 'VNP02IMG',
    #     'enhancers': 'default'}],
    # 'temporal_interp': 'linear',
    # 'spatial_interp': 'nearest'}

    adjust_logger(True, folder, "INFO")

    ds = pywapor.pre_se_root.main(folder, latlim, lonlim, timelim, sources)

    ds_out = pywapor.se_root.main(ds, se_root_version = "v3")