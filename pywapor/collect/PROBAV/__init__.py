# -*- coding: utf-8 -*-
"""
Authors: Laust Færch
Module: Collect/PROBAV

Description:
This module downloads Proba-V NDVI data from
https://www.vito-eodata.be/PDF/datapool/.
You need to register as a user (which is free) to use this module

"""
from inspect import getsourcefile
import sys
from pathlib import Path

probav_init_path = getsourcefile(lambda:0)
probav_module_path = Path(probav_init_path).parent

vito_download_path = Path.joinpath(probav_module_path, "vito-download")
itsybitsy_path = Path.joinpath(probav_module_path, "itsybitsy")

# Add paths for Gitmodules to sys.path
sys.path.append(vito_download_path.as_posix())
sys.path.append(itsybitsy_path.as_posix())

# print(sys.path)

from .PROBAV_S5 import main as PROBAV_S5

__all__ = ["PROBAV_S5"]

__version__ = '0.1'
