# -*- coding: utf-8 -*-
"""
Generates input data for `pywapor.se_root`.
"""
from pywapor.collect import downloader
from pywapor.general import compositer
from pywapor.general import aligner

def main(folder, latlim, lonlim, timelim, sources, example_source, bin_length):

    bins = compositer.time_bins(timelim, bin_length)
    dss = downloader.collect_sources(folder, sources, latlim, lonlim, [bins[0], bins[-1]])
    ds = aligner.main(dss, sources, example_source, folder)

    return ds

if __name__ == "__main__":

    import datetime

    sources = {

        "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
        "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],
        "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        "u2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        "v2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        "qv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        "wv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        "p_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        "p_air_0":      [("MERRA2", "M2I1NXASM.5.12.4")],

    }

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 6, 25), datetime.date(2020, 7, 30)]
    example_source = ("MODIS", "MOD13Q1.061")
    bin_length = 4

    ds = main(folder, latlim, lonlim, timelim, sources, example_source, bin_length)