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

    total = 0
    for sd, ed, periods in zip(Startdate, Enddate, Periods):
        no = len(pd.date_range(sd, ed, freq = "D")) * len(periods)
        total += no

    super_total = total * no_vars

    if len(Vars) == 1:
        log.info(f"--> Downloading GEOS5 (3-hourly), {Vars[0]}.")
    else:
        log.info(f"--> Downloading GEOS5 (3-hourly), {Vars}.")

    if Waitbar:
        # waitbar = tqdm.tqdm(total = no_vars * no_dates, position = 0, unit = "tiles")

        waitbar = tqdm.tqdm(desc= f"Tile: 0 / {super_total}",
                            position = 0,
                            # total=total_size,
                            unit='Bytes',
                            unit_scale=True,)

    all_files = list()

    for Var in Vars:

        for sd, ed, periods in zip(Startdate, Enddate, Periods):

            for Period in periods:

                # Download data
                files = DownloadData(Dir, Var, sd, ed, latlim, lonlim, "three_hourly", Period, waitbar)
                all_files = all_files + files


    if Waitbar:
        waitbar.close()

    all_files = np.unique(all_files).tolist()

    return all_files

if __name__ == '__main__':
    
    Dir = r"/Users/hmcoerver/pywapor_notebooks/RAW"
    latlim = [-40.0, 20.0]
    lonlim = [80.0, 100.0]
    Startdate = ["2021-07-01"]
    Enddate = ["2021-07-03"]
    Periods = [[1,2,3]]

    Vars = ['u2m']

    main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Periods, Waitbar = True)

