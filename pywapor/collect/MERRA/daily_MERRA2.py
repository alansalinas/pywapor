# -*- coding: utf-8 -*-
import sys
from pywapor.collect.MERRA.DataAccess import DownloadData
from pywapor.collect.MERRA.DataAccess import VariablesInfo
import pywapor

def main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Waitbar = 1, data_type = ["mean"]):
    """
    This function downloads MERRA daily data for a given variable, time
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
    all_files = dict()

    username, password = pywapor.collect.get_pw_un.get("NASA")

    for Var in Vars:

        if Waitbar == 1:
            print('\nDownloading daily MERRA %s data for the period %s till %s' %(Var, Startdate, Enddate))

        # Download data
        all_files = {**all_files, **DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "daily_MERRA2", '', username, password, Waitbar = 1, data_type = data_type)}

        TimeStep = "daily_MERRA2"
        VarInfo = VariablesInfo(TimeStep)
        Parameter = VarInfo.names[Var]

    return all_files

if __name__ == '__main__':
    main(sys.argv)
