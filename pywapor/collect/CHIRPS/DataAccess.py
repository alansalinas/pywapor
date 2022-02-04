# -*- coding: utf-8 -*-
"""
Authors: Tim Hessels
Module: Collect/CHIRPS
"""

# import general python modules
import os
from osgeo import gdal
import numpy as np
import pandas as pd
from ftplib import FTP
import datetime
import tqdm

def DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar):
    """
    This function downloads CHIRPS daily or monthly data

    Keyword arguments:
    Dir -- 'C:/file/to/path/'
    Startdate -- 'yyyy-mm-dd'
    Enddate -- 'yyyy-mm-dd'
    latlim -- [ymin, ymax] (values must be between -50 and 50)
    lonlim -- [xmin, xmax] (values must be between -180 and 180)
    Waitbar -- 1 (Default) will print a waitbar
    cores -- The number of cores used to run the routine. It can be 'False'
             to avoid using parallel computing routines.
    TimeCase -- String equal to 'daily' or 'monthly'
    """
	
    # Define timestep for the timedates
    output_folder = os.path.join(Dir, "CHIRPS", "Precipitation")

    # make directory if it not exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

	# check time variables
    if not Startdate:
        Startdate = pd.Timestamp('1981-01-01')
    if not Enddate:
        Enddate = pd.Timestamp('Now')

    # Create days
    Dates = pd.date_range(Startdate, Enddate, freq='D')

    if Waitbar:
        waitbar = tqdm.tqdm(desc= f"Tile: 0 / {int(len(Dates))}",
                            position = 0,
                            # total=total_size,
                            unit='Bytes',
                            unit_scale=True,)
    else:
        waitbar = None

    # Check space variables
    if latlim[0] < -50 or latlim[1] > 50:
        print('Latitude above 50N or below 50S is not possible.'
               ' Value set to maximum')
        latlim[0] = np.max(latlim[0], -50)
        latlim[1] = np.min(lonlim[1], 50)
    if lonlim[0] < -180 or lonlim[1] > 180:
        print('Longitude must be between 180E and 180W.'
               ' Now value is set to maximum')
        lonlim[0] = np.max(latlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)

    # Define IDs
    yID = 2000 - np.int16(np.array([np.ceil((latlim[1] + 50)*20),
                                    np.floor((latlim[0] + 50)*20)]))
    xID = np.int16(np.array([np.floor((lonlim[0] + 180)*20),
                             np.ceil((lonlim[1] + 180)*20)]))

    # Pass variables to parallel function and run
    args = [output_folder, xID, yID, lonlim, latlim]
    for Date in Dates:
        RetrieveData(Date, args, waitbar)

    results = True
    
    # remove raw files
    import glob
    os.chdir(output_folder)
    try:
        files_raw = glob.glob("chirps-v2.0*.tif")
        for file_raw in files_raw:
            outfilename = os.path.join(output_folder, file_raw)
            # delete old tif file
            os.remove(outfilename)
    except:
        pass        
        
    return results

def RetrieveData(Date, args, waitbar):
    """
    This function retrieves CHIRPS data for a given date from the
    ftp://chg-ftpout.geog.ucsb.edu server.

    Keyword arguments:
    Date -- 'yyyy-mm-dd'
    args -- A list of parameters defined in the DownloadData function.
    """
    
    # WA+ modules
    import pywapor.general.processing_functions as PF  

    # Argument
    [output_folder, xID, yID, lonlim, latlim] = args

    # Define output
    DirFileEnd = os.path.join(output_folder,'P_CHIRPS.v2.0_mm-day-1_daily_%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))

    if not os.path.exists(DirFileEnd):

        # open ftp server
        ftp = FTP("ftp.chc.ucsb.edu", "", "")
        ftp.login()
    
    	# Define FTP path to directory
        pathFTP = 'pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/p05/%s/' %Date.strftime('%Y')
    
        # find the document name in this directory
        ftp.cwd(pathFTP)
        listing = []
    
    	# read all the file names in the directory
        ftp.retrlines("LIST", listing.append)
    
    	# create all the input name (filename) and output (outfilename, filetif, DiFileEnd) names
        filename = 'chirps-v2.0.%s.%02s.%02s.tif.gz' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d'))
        outfilename = os.path.join(output_folder,'chirps-v2.0.%s.%02s.%02s.tif' %(Date.strftime('%Y'), Date.strftime('%m'), Date.strftime('%d')))

        # download the global rainfall file

        local_filename = os.path.join(output_folder, filename)
        total_size = ftp.size(filename)
        if not isinstance(waitbar, type(None)):
            waitbar_i = int(waitbar.desc.split(" ")[1])
            waitbar_desc = str(waitbar.desc)
            waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /"))
            waitbar.reset(total = total_size)
        with open(local_filename, "wb") as lf:
            def callback(data):
                lf.write(data)
                if not isinstance(waitbar, type(None)):
                    waitbar.update(len(data))
            ftp.retrbinary("RETR " + filename, callback, 8192) 
        PF.Extract_Data_gz(local_filename, outfilename)

        # open tiff file
        dest = gdal.Open(outfilename)
        dataset = dest.GetRasterBand(1).ReadAsArray()

        # clip dataset to the given extent
        data = dataset[yID[0]:yID[1], xID[0]:xID[1]]
        data[data < 0] = -9999

        # save dataset as geotiff file
        geo = [lonlim[0], 0.05, 0, latlim[1], 0, -0.05]
        
        PF.Save_as_tiff(name=DirFileEnd, data=data, geo=geo, projection="WGS84")
    else:
        if not isinstance(waitbar, type(None)):
            waitbar_i = int(waitbar.desc.split(" ")[1])
            waitbar_desc = str(waitbar.desc)
            waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /"))

    return True

if __name__ == "__main__":

    import pywapor

    Dir =r"/Volumes/Data/FAO/temp_test"
    Startdate = "2021-06-01"
    Enddate = "2021-07-01"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]

    DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar = 1)
