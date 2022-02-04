# pyWAPOR
![downloads](https://img.shields.io/pypi/dw/pywapor) [![version](https://img.shields.io/pypi/v/pywapor)](https://pypi.org/project/pywapor/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/pywapor_101.ipynb)

This repository contains a Python implementation of the algorithm used to generate the [WaPOR](http://www.fao.org/in-action/remote-sensing-for-water-productivity/en/) [datasets](https://wapor.apps.fao.org/home/WAPOR_2/1). It can be used to calculate evaporation, transpiration and biomass production maps.

## Installation

Its recommended to install in a clean [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install GDAL before installing pywapor.

```bash
conda create -n my_pywapor_env python=3.8 pip gdal
conda activate my_pywapor_env
```

On machines running Windows also run the following.

```bash
conda install -c conda-forge fiona rasterio pycurl
```

Then use the package manager [pip](https://pip.pypa.io/en/stable/) to install pywapor.

```bash
pip install pywapor
```

## Usage

To run the model for one dekad (from 2021-07-01 to 2021-07-11 in this case) for the Fayoum irrigation scheme in Egypt (but feel free to change the [boundingbox](http://bboxfinder.com) defined by `latlim` and `lonlim`) using mainly MODIS data, run the following code. 

```python
import pywapor

# User inputs.
startdate = "2021-07-01"
enddate = "2021-07-11"
latlim = [28.9, 29.7]
lonlim = [30.2, 31.2]
project_folder = r"/my_first_ETLook_run/"

# Download input data.
ds_in, fh_in = pywapor.pre_et_look.main(project_folder, startdate, enddate, latlim, lonlim)

# Run the model for one dekad starting on 'startdate'.
ds_out = pywapor.et_look.main(ds_in)
```

Check out one of the Colab Notebooks below to learn more!

### Notebooks
|  | Name | Duration* |
| ------ | ------ | ------ |
| 1. | [Introduction](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/1_introduction.ipynb) | 10 + 30 |
| 2. | [Composites](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/2_composites.ipynb) | 10 + 10 |
| 3. | [Levels](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/3_levels.ipynb) | 10 + 120 |
| 4. | [Sideloading](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/4_sideloading.ipynb) | 10 + 5 |
| 5. | [Enhancers](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/5_enhancers.ipynb) | 10 + 5 |
| 6. | [pyWaPOR vs. WaPOR](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/6_wapor_vs_pywapor.ipynb) | 10 + 30 |
| 7. | Soil Saturation | > Coming Soon < |

\* Estimation of the time required in minutes, as in "active" + "download time". 

## Documentation
### WaPOR v2
➡ [WaPOR-ETLook Data Manual](https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_ET_v2_data_manual_finaldraft_v2.2.pdf)

➡ [WaPOR-Biomass Data Manual](https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_NPP_v2_data_manual_finaldraft_v2.2.pdf)

### WaPOR v1
➡ [WaPOR-ETLook Data Manual](https://bitbucket.org/cioapps/wapor-et-look/raw/9ec88e56769f49722c2d1165bb34547f5842b811/Docs/WaPOR_ET_data_manual_finaldraft-v1.2-for-distribution.pdf)

## Acknowledgments
The methodology for WaPOR was developed by the FRAME1 consortium, consisting of eLEAF, VITO, ITC, University of Twente and Waterwatch foundation, commissioned by and in partnership with the Land and Water Division of FAO. 

This repository contains, among others, contributions from Bert Coerver (FAO), Tim Hessels (WaterSat) and, in the framework of the ESA-funded ET4FAO project, from Radoslaw Guzinski (DHI-GRAS), Hector Nieto (Complutig) and Laust Faerch (DHI-GRAS).

## Contact
For questions, requests or issues with this repository, please contact Bert Coerver at [bert.coerver@fao.org](mailto:bert.coerver@fao.org) or the WaPOR team at [wapor@fao.org](mailto:wapor@fao.org).

## Release Notes

#### 2.4.0 (2022-02-03)

* Easily apply your own functions to data, i.e. use your own custom filters, gap-fillers etc.
* Side-load your own data, i.e. easily incorporate you own datasets.
* Added functions to process Landsat Level-2 to “ndvi”, “lst” and “r0”.
* Data is now stored and processed as netCDF (using xarray and dask).
* et_look() and se_root() now calculate in chunks, instead of using a for-loop.
* Some previously constant parameters now have spatial variability.
* Improved logging.
* Download functions now show progress and download-speed.
* MODIS datasets switched from v6.0 (decommissioned soon) to v6.1.
* The lapse-rate correction to temperature data is now more accurate and faster.
* VITO and WAPOR passwords are now checked when entered.
* Other bug-fixes and performance improvements.

#### 2.3.0 (2021-11-19)

* Automatically create input composites before running ETLook.
* Choose composite lengths in number of days or dekads.
* Option to choose which products to use per variable.
* Calculate soil saturation separate from ETLook.
* PROBA-V support for NDVI and Albedo inputs.
* Define diagnostics pixels, for which extra outputs are created (e.g. charts, maps etc.).
* Bug-fixes and performance improvements.

## License
[APACHE](https://bitbucket.org/cioapps/wapor-et-look/src/dev/LICENSE)