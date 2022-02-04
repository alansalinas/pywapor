import sys
from pywapor.collect.MYD13.DataAccess import DownloadData
from datetime import date, timedelta
import numpy as np
from datetime import datetime as dat
import glob
import os
import pywapor
from pywapor.general.logger import log

def main(Dir, latlim, lonlim, Startdate, Enddate, Waitbar = 1, 
        hdf_library = None, remove_hdf = 1, buffer_dates = False):
    """
    This function downloads MYD13 16-daily data for the specified time
    interval, and spatial extent.

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd' or datetime.date
    Enddate -- 'yyyy-mm-dd' or datetime.date
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    username -- "" string giving the username of your NASA account (https://urs.earthdata.nasa.gov/)
    password -- "" string giving the password of your NASA account    
    Waitbar -- 1 (Default) will print a waitbar
    hdf_library -- string, if all the hdf files are already stored on computer
                    define directory to the data here
    remove_hdf -- 1 (Default), if 1 remove all the downloaded hdf files in the end    
    """
    if isinstance(Startdate, date):
        Startdate = Startdate.strftime("%Y-%m-%d")
    
    if isinstance(Enddate, date):
        Enddate = Enddate.strftime("%Y-%m-%d")

    username, password = pywapor.collect.get_pw_un.get("NASA")

    log.info(f"--> Downloading MYD13.")
    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, username, 
                password, Waitbar, hdf_library, remove_hdf, buffer_dates = buffer_dates)

    sdate = dat.strptime(Startdate, "%Y-%m-%d")
    edate = dat.strptime(Enddate, "%Y-%m-%d")
    all_files = glob.glob(os.path.join(Dir, "MODIS", "MYD13", "*.tif"))
    start_dates = np.array([dat.strptime(os.path.split(x)[-1], "NDVI_MYD13Q1_-_16-daily_%Y.%m.%d.tif") for x in all_files])
    end_dates = np.array([x + timedelta(days = 16) for x in start_dates])
    check = np.all([start_dates <= edate, end_dates >= sdate], axis = 0)
    all_files = np.array(all_files)[check].tolist()

    return all_files


# if __name__ == '__main__':
#     main(sys.argv)