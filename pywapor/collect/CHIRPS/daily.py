# -*- coding: utf-8 -*-
import sys
from pywapor.collect.CHIRPS.DataAccess import DownloadData
import glob
import os
from pywapor.general.logger import log

def main(Dir, latlim, lonlim, Startdate, Enddate, Waitbar = 1):
	"""
	This function downloads daily CHIRPS data

	Keyword arguments:
	Dir -- 'C:/file/to/path/'
	Startdate -- 'yyyy-mm-dd'
	Enddate -- 'yyyy-mm-dd'
	latlim -- [ymin, ymax] (values must be between -50 and 50)
	lonlim -- [xmin, xmax] (values must be between -180 and 180)
	Waitbar -- 1 (Default) will print a waitbar
	"""

	log.info(f"--> Downloading CHIRPS.")
	# Download data
	DownloadData(Dir, Startdate, Enddate, latlim, lonlim, Waitbar)

	return glob.glob(os.path.join(Dir, "CHIRPS", "Precipitation", "*.tif"))

if __name__ == '__main__':
    main(sys.argv)
