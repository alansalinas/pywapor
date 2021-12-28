# -*- coding: utf-8 -*-
import sys
from pywapor.collect.GEOS5.DataAccess import DownloadData
import pandas as pd
import tqdm
from pywapor.general.logger import log
import datetime
import numpy as np

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Periods = [1,2,3,4,5,6,7,8], Waitbar = True):
    """
    This function downloads GEOS inst data for a given variable, time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Vars -- ['t2m', 'v2m'] 
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
	Periods -- [1,2,3,4,5,6,7,8] Period that needs to be downloaded. 1 period is 3 hour starting from noon
    Waitbar -- 1 (Default) Will print a waitbar
    """
    no_vars = len(Vars)

    if isinstance(Waitbar, tqdm.tqdm):
        waitbar = Waitbar
    elif Waitbar:
        log.info("--> Downloading GEOS5 (hourly).")
        waitbar = tqdm.tqdm(total = no_vars * sum([len(x) for x in Periods]))
    else:
        log.info("--> Downloading GEOS5 (hourly).")
        waitbar = False

    all_files = list()

    for Var in Vars:

        if Waitbar:
            waitbar.set_description_str(Var)

        for sd, ed, periods in zip(Startdate, Enddate, Periods):

            for Period in periods:
            
                # if Waitbar == 1:
                #     print('\nDownloading 3-hourly GEOS %s data for the period %s till %s, Period = %s' %(Var, Startdate, Enddate, Period))

                # Download data
                files = DownloadData(Dir, Var, sd, ed, latlim, lonlim, "three_hourly", Period, waitbar)
                all_files = all_files + files


    if Waitbar:
        waitbar.close()

    all_files = np.unique(all_files).tolist()

    return all_files

# if __name__ == '__main__':
#     main(sys.argv)
