# -*- coding: utf-8 -*-
"""
Authors: Laust Færch
Module: Collect/PROBAV

Description:
This module downloads Proba-V NDVI data from
https://www.vito-eodata.be/PDF/datapool/.
You need to register as a user (which is free) to use this module

"""
import os
from .PROBAV_S5 import main as PROBAV_S5
import inspect
from pathlib import Path
# import pywapor.collect.PROBAV as probav

probav_module_path = Path(inspect.getfile(PROBAV_S5)).parent
sub_module_path = Path.joinpath(probav_module_path, "vito-download")
current_workdir = os.getcwd()
os.chdir(sub_module_path)
import vito_download as vito
os.chdir(current_workdir)

print(vito)

__all__ = ["PROBAV_S5"]

__version__ = '0.1'
