# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/Globcover
"""
from pywapor.collect.Globcover.DataAccess import DownloadData
import sys
import os

def main(Dir, latlim, lonlim, Waitbar = 1):
    """
    Downloads Globcover data from http://due.esrin.esa.int/page_globcover.php

    The following keyword arguments are needed:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- '1' if you want a waitbar (Default = 1)
    """
    
    # Create Waitbar
    if Waitbar == 1:
        print('\nDownload Globcover landcover map')
        import pywapor.general.waitbar_console as WaitbarConsole
        total_amount = 1
        amount = 0
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    # Download and process the data
    DownloadData(Dir, latlim, lonlim)

    if Waitbar == 1:
        amount = 1
        WaitbarConsole.printWaitBar(amount, total_amount, prefix = 'Progress:', suffix = 'Complete', length = 50)

    outpath = os.path.join(Dir, "GLOBCOVER", "Landuse", "LC_GLOBCOVER_V2.3.tif")

    return outpath

if __name__ == '__main__':
    main(sys.argv)
