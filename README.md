# pyWAPOR
![downloads](https://img.shields.io/pypi/dw/pywapor)
[![version](https://img.shields.io/pypi/v/pywapor)](https://pypi.org/project/pywapor/)

This repository contains a Python implementation of the algorithm used to generate the [WaPOR](http://www.fao.org/in-action/remote-sensing-for-water-productivity/en/) [datasets](https://wapor.apps.fao.org/home/WAPOR_2/1). It can be used to calculate evaporation, transpiration and biomass production maps based on MODIS or Landsat data.

## Installation

Its recommended to install in a clean [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install GDAL before installing pywapor.

```bash
conda create -n my_pywapor_env python=3 pip
conda activate my_pywapor_env
conda install -c conda-forge gdal=3.1.4
```

Then use the package manager [pip](https://pip.pypa.io/en/stable/) to install pywapor.

```bash
pip install pywapor
```

## Usage

```python
import os
import pywapor

# User inputs.
startdate = "2019-07-07"
enddate = "2019-07-07"
latlim = [29.0, 29.6]
lonlim = [30.3, 31.1]
project_folder = r"/my_first_ETLook_run/"

# Download input data.
pywapor.pre_et_look.main(project_folder, startdate, enddate, latlim, lonlim)

# Run the model.
ETLook_input_folder = os.path.join(project_folder, "ETLook_input_MODIS")
ETLook_output_folder = os.path.join(project_folder, "ETLook_output_MODIS")

pywapor.et_look_code.main(ETLook_input_folder, ETLook_output_folder, startdate)
```

See the examples folder for more examples or check out the [Colab Notebook](https://colab.research.google.com/drive/1YEsCN6GnMGvOzXT4YaJ_jeu58mGIhhMq?usp=sharing).

## Documentation
### WaPOR v2
➡ [WaPOR-ETLook Data Manual](https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_ET_v2_data_manual_finaldraft_v2.2.pdf)

➡ [WaPOR-Biomass Data Manual](https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_NPP_v2_data_manual_finaldraft_v2.2.pdf)

### WaPOR v1
➡ [WaPOR-ETLook Data Manual](https://bitbucket.org/cioapps/wapor-et-look/raw/9ec88e56769f49722c2d1165bb34547f5842b811/Docs/WaPOR_ET_data_manual_finaldraft-v1.2-for-distribution.pdf)

## Acknowledgments
The methodology for WaPOR was developed by the FRAME1 consortium, consisting of eLEAF, VITO, ITC, University of Twente and Waterwatch foundation, commissioned by and in partnership with the Land and Water Division of FAO.

## License
[APACHE](https://bitbucket.org/cioapps/wapor-et-look/src/dev/LICENSE)
