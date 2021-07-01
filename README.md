# pyWAPOR

This repository contains a Python implementation of the algorithm used to generate the [WaPOR](http://www.fao.org/in-action/remote-sensing-for-water-productivity/en/) [datasets](https://wapor.apps.fao.org/home/WAPOR_2/1). It can be used to calculate evaporation, transpiration and biomass production maps based on MODIS or Landsat data.

## Installation

Its recommended to install in a clean [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install GDAL before installing pywapor.

```bash
conda create -n my_pywapor_env python=3 pip
conda activate my_pywapor_env
conda install -c conda-forge gdal
```

Then use the package manager [pip](https://pip.pypa.io/en/stable/) to install pywapor.

```bash
pip install pywapor
```

## Usage

```python
import os
import pandas as pd
import pywapor

# User Inputs.
startdate = "2019-07-06"
enddate = "2019-07-06"
latlim = [28.5, 31.9]
lonlim = [29.2, 32.5]
output_folder = r"/path/to/local/folder"
etlook_version = "v2"
##############

dates = pd.date_range(startdate, enddate, freq = "D")
raw_folder = os.path.join(output_folder, "RAW")
model_input = os.path.join(output_folder, "ETLook_input")
model_output = os.path.join(output_folder, "ETLook_output")

# Download input data.
pywapor.pre_et_look.main(
                         output_folder, 
                         startdate, 
                         enddate, 
                         latlim, 
                         lonlim,
                         RAW_folder= raw_folder,
                       )

# Run the model.
for date in dates:
    pywapor.et_look_code.main(model_input, model_output,
                              dates, ETLook_version = etlook_version)
```

See the examples folder for more examples or check out the [Colab Notebook]().

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
