# -*- coding: utf-8 -*-
from pywapor.collect.GEOS5.DataAccess import DownloadData
import tqdm
import pandas as pd
from pywapor.general.logger import log

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

    if len(Vars) == 1:
        log.info(f"--> Downloading GEOS5 (daily), {Vars[0]}.")
    else:
        log.info(f"--> Downloading GEOS5 (daily), {Vars}.")

    if Waitbar:
        # waitbar = tqdm.tqdm(total = no_vars * no_dates, position = 0, unit = "tiles")

        waitbar = tqdm.tqdm(desc= f"Tile: 0 / {no_vars * no_dates}",
                            position = 0,
                            # total=total_size,
                            unit='Bytes',
                            unit_scale=True,)

    else:
        waitbar = False

    all_files = list()
    for Var in Vars:

        if Var == "t2m":
            filter = "t2m_"
        elif Var == "t2m-max":
            Var = "t2m"
            filter = "t2m-max_"
        elif Var == "t2m-min":
            Var = "t2m"
            filter = "t2m-min_"
        else:
            filter = None

        # if Waitbar:
        #     waitbar.set_description_str("{:<11}".format(Var))

        # Download data
        new_files = DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, "daily", Period = '', Waitbar = waitbar)
        
        if not isinstance(filter, type(None)):
            new_files = [x for x in new_files if filter in x]
        
        all_files = all_files + new_files

    if Waitbar:
        waitbar.close()

    return all_files

if __name__ == '__main__':
    
    Dir = r"/Users/hmcoerver/pywapor_notebooks/RAW"
    latlim = [-40.0, 20.0]
    lonlim = [80.0, 100.0]
    Startdate = "2021-07-01"
    Enddate = "2021-07-03"

    Vars = ['u2m']

    main(Dir, latlim, lonlim, Startdate, Enddate, Vars, Waitbar = True)
