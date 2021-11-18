# -*- coding: utf-8 -*-
import sys
from pywapor.collect.GEOS.DataAccess import DownloadData
import tqdm
import pandas as pd

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Waitbar = True):
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
    no_dates = len(pd.date_range(Startdate, Enddate, freq="D"))
    no_vars = len(Vars)

    if Waitbar:
        waitbar = tqdm.tqdm(total = no_vars * no_dates, delay = 30, position = 0)
    else:
        waitbar = False

    all_files = dict()
    for Var in Vars:

        if Waitbar:
            waitbar.set_description_str(Var)

        # Download data
        all_files = {**all_files, **DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "daily", Period = '', Waitbar = waitbar)}


    return all_files

if __name__ == '__main__':
    main(sys.argv)
