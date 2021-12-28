# -*- coding: utf-8 -*-
import sys
from pywapor.collect.MERRA2.DataAccess import DownloadData
import pywapor
from pywapor.general.logger import log
import numpy as np

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Periods = list(range(1, 25)), Waitbar = False, verbose = True):
    """
    This function downloads MERRA inst data for a given variable, time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Vars -- ['t2m', 'v2m'] 
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
	Periods -- [1,2,3,4,5,6,7,8,23,24] Period that needs to be downloaded. 1 period is 1 hour starting from noon
    Waitbar -- 1 (Default) Will print a waitbar
    """

    username, password = pywapor.collect.get_pw_un.get("NASA")

    all_files = list()

    if not isinstance(Startdate, list):
        Startdate = [Startdate]
    if not isinstance(Enddate, list):
        Startdate = [Enddate]

    if not isinstance(Periods[0], list):
        Periods = [Periods for _ in Startdate]
    
    if not verbose:
        log.info("--> Downloading MERRA2 (hourly).")
    
    for Var in Vars:

        for sd, ed, periods in zip(Startdate, Enddate, Periods):

            for period in periods:

                # Download data
                files = DownloadData(Dir, Var, sd, ed, latlim, lonlim, "hourly_MERRA2", period, username, password, Waitbar, data_type = ["mean"])
                all_files = all_files + files

    all_files = np.unique(all_files).tolist()

    return all_files

if __name__ == '__main__':
    main(sys.argv)

