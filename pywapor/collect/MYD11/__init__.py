# -*- coding: utf-8 -*-
"""

Description:
This module downloads MYD11 LST data from
http://e4ftl01.cr.usgs.gov/. Use the MYD11.LST function to
download and create daily LST images in Gtiff format.
The data is available between 2000-02-18 till present.

Examples:
from pyWAPOR.Collect import MYD11
MYD11.LST(Dir='C:/Temp3/', Startdate='2003-12-01', Enddate='2003-12-30',
           latlim=[41, 45], lonlim=[-8, -5])
"""

from pywapor.collect.MYD11.LST import main as LST

__all__ = ['LST']

__version__ = '0.1'
