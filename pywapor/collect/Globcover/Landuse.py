from pywapor.collect.Globcover.DataAccess import DownloadData
import sys
import os
from pywapor.general.logger import log

def main(Dir, latlim, lonlim, Waitbar = 1, **kwargs):
    """
    Downloads Globcover data from http://due.esrin.esa.int/page_globcover.php

    The following keyword arguments are needed:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]

    """
    log.info("--> Downloading GLOBCOVER.")

    # Download and process the data
    outpath = DownloadData(Dir, latlim, lonlim)

    return outpath

if __name__ == '__main__':
    main(sys.argv)
