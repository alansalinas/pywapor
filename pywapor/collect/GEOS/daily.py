# -*- coding: utf-8 -*-
import sys
from pywapor.collect.GEOS.DataAccess import DownloadData
import os
import glob

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Waitbar = 1):
    """
    This function downloads GEOS daily data for a given variable, time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Vars -- ['t2m', 'v2m'] 
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    Waitbar -- 1 (Default) Will print a waitbar
    """
    for Var in Vars:

        if Waitbar == 1:
            print('\nDownloading daily GEOS %s data for the period %s till %s' %(Var, Startdate, Enddate))

        # Download data
        DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "daily", Period = '', Waitbar = 1)

    output_folder = os.path.join(Dir, "GEOS5", "**", "**", "*.tif")
    return glob.glob(output_folder)

if __name__ == '__main__':
    main(sys.argv)
