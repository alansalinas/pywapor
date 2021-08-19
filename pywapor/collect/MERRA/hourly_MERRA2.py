# -*- coding: utf-8 -*-
import sys
from pywapor.collect.MERRA.DataAccess import DownloadData
import pywapor

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Periods = list(range(1, 25)), Waitbar = 1):
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

    for Var in Vars:

        for Period in Periods:
		
            if Waitbar == 1:
                print('\nDownloading hourly MERRA %s data for the period %s till %s, Period = %s' %(Var, Startdate, Enddate, Period))

            # Download data
            DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "hourly_MERRA2", Period, username, password, Waitbar, data_type = ["mean"])

if __name__ == '__main__':
    main(sys.argv)

