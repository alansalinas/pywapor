# -*- coding: utf-8 -*-
"""

This module downloads daily CHIRPS 2.0 data from
ftp://chg-ftpout.geog.ucsb.edu server. The CHIRPS data is available since 1981-01-01 till the present.
The datasets will be stored in the user defined outputfolder in GEOTIFF format.


Examples:
from pyWAPOR.Collect import CHIRPS
CHIRPS.daily(Dir='C:/Temp/', Startdate='1999-02-01', Enddate='1999-02-03',
             latlim=[-10, 30], lonlim=[-20, 120])
"""

from pywapor.collect.CHIRPS.PRECIPITATION import main as PRECIPITATION

__all__ = ['PRECIPITATION']

__version__ = '0.1'
