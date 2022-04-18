# -*- coding: utf-8 -*-
"""
Generates input data for `pywapor.et_look`.
"""
from pywapor.collect import downloader
from pywapor.general.logger import log, adjust_logger
from pywapor.general import compositer

def main(folder, latlim, lonlim, timelim, sources, example_source, 
                    bin_length, diagnostics = None):
    
    bins = compositer.time_bins(timelim, bin_length)
    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])
    ds = compositer.main(dss, sources, example_source, bins, folder)
    
    if diagnostics:
        _ = compositer.main(dss, sources, example_source, bins, folder, diagnostics = diagnostics)
    
    return ds

if __name__ == "__main__":

    import datetime

    sources = {
        "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
        "r0":           [("MODIS", "MCD43A3.061")],
        "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
        "z":            [("SRTM", "30M")],
        "p":            [("CHIRPS", "P05")],
        "ra":           [("MERRA2", "M2T1NXRAD.5.12.4")],
        "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        # "t_air_max":    [("MERRA2", "M2I1NXASM.5.12.4")],
    }

    diagnostics = { # label          # lat      # lon
                    "water":	    (29.44977,	30.58215),
                    "desert":	    (29.12343,	30.51222),
                    "agriculture":	(29.32301,	30.77599),
                    "urban":	    (29.30962,	30.84109),
                    }

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
    example_source = ("MODIS", "MOD13Q1.061")
    bin_length = 4

    ds = main(folder, latlim, lonlim, timelim, sources, example_source, 
                bin_length, diagnostics = diagnostics)


