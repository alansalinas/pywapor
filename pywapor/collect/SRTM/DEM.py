
import os
from pywapor.collect.SRTM.DataAccess import DownloadData
import sys


def main(Dir, latlim, lonlim, **kwargs):
    """
    Downloads HydroSHED data from http://srtm.csi.cgiar.org/download

    this data includes a Digital Elevation Model (DEM)
    The spatial resolution is 90m (3s)

    The following keyword arguments are needed:
    Dir -- 'C:/file/to/path/'
    latlim -- [ymin, ymax]
    lonlim -- [xmin, xmax]
    """

    # Create directory if not exists for the output
    output_folder = os.path.join(Dir, 'SRTM', 'DEM')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define the output map and create this if not exists
    nameEnd = os.path.join(output_folder, 'DEM_SRTM_m_3s.tif')

    if not os.path.exists(nameEnd):

        # Download and process the data
        DownloadData(output_folder, latlim, lonlim)

    return nameEnd

if __name__ == '__main__':
    main(sys.argv)
