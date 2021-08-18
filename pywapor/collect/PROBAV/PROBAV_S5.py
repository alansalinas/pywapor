import sys
from pathlib import Path
from pywapor.collect.PROBAV.DataAccess import download_data


def main(download_dir, start_date, end_date, latitude_extent, longitude_extent, username,
         password, buffer_dates = True):
    """

    """

    download_data(download_dir/Path("PROBAV"), start_date, end_date, latitude_extent,
                  longitude_extent, username, password, buffer_dates = buffer_dates)


if __name__ == '__main__':
    main(sys.argv)
