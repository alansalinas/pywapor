# -*- coding: utf-8 -*-
import sys
from pywapor.collect.GEOS.DataAccess import DownloadData
import pandas as pd
import tqdm

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
    no_periods = len(Periods)
    no_dates = len(pd.date_range(Startdate, Enddate, freq="D"))

    if isinstance(Waitbar, tqdm.tqdm):
        waitbar = Waitbar
    elif Waitbar:
        waitbar = tqdm.tqdm(total = no_vars * no_dates * no_periods, delay = 10)
    else:
        waitbar = False

    for Var in Vars:

        if Waitbar:
            waitbar.set_description_str(Var)

        for Period in Periods:
		
            # if Waitbar == 1:
            #     print('\nDownloading 3-hourly GEOS %s data for the period %s till %s, Period = %s' %(Var, Startdate, Enddate, Period))

            # Download data
            DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "three_hourly", Period, waitbar)

if __name__ == '__main__':
    main(sys.argv)
