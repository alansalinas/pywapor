# pyWAPOR
![downloads](https://img.shields.io/pypi/dw/pywapor) [![version](https://img.shields.io/pypi/v/pywapor)](https://pypi.org/project/pywapor/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/pywapor_101.ipynb)


This repository contains a Python implementation of the algorithm used to generate the [WaPOR](http://www.fao.org/in-action/remote-sensing-for-water-productivity/en/) [datasets](https://wapor.apps.fao.org/home/WAPOR_2/1). It can be used to calculate evaporation, transpiration and biomass production maps.

## Installation

Its recommended to install in a clean [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install GDAL before installing pywapor.

```bash
conda create -n my_pywapor_env python=3 pip
conda activate my_pywapor_env
conda install -c conda-forge gdal=3.1.4
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

To run the model for one dekad (from 2021-07-01 to 2021-07-10 in this case) for the Fayoum irrigation scheme in Egypt (but feel free to change the [boundingbox](http://bboxfinder.com) defined by `latlim` and `lonlim`) using mainly MODIS data, run the following code. 

```python
import pywapor

# User inputs.
startdate = "2021-07-01"
enddate = "2021-07-10"
latlim = [28.9, 29.7]
lonlim = [30.2, 31.2]
project_folder = r"/my_first_ETLook_run/"

# Download input data.
pywapor.pre_et_look.main(project_folder, startdate, enddate, latlim, lonlim)

# Run the model for one day.
pywapor.et_look.main(project_folder, startdate)
```

Check out one of the Colab Notebooks below to learn more!

### Notebooks
|  | Name | Duration* |
| ------ | ------ | ------ |
| 1. | [Introduction](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/pywapor_101.ipynb) | 10 + 30 |
| 2. | [Levels](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/levels.ipynb) | 10 + 120 |
| 3. | [Composites](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/composites.ipynb) | 10 + 10 |
| 4. | [pyWaPOR vs. WaPOR](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/wapor_vs_pywapor.ipynb) | 10 + 30 |
| 5. | Soil Saturation | > Coming Soon < |

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

## Data Sources

To run the ETLook model, two types of spatial variables are required, temporal and static data. **Each of these variables can be collected from whichever source you wish to use**, as long as you make sure the units are correct, the data is stored as a GeoTIFF (1 band per file, 1 file for each variable and date), the files all have the same no-data-value and they all have the same projection and resolution.

**For your convenience, the pyWAPOR package has a function that can collect all this data from selected sources** and make sure the data is stored in the correct format and folder structure.

#### Temporal ET_Look Data (composites)
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Normalized Difference Vegetation Index (NDVI) | - | MOD13, MYD13, PROBA-V|
| Albedo | - | MCD43, PROBA-V|
| Precipitation | mm/day | CHIRPS |
| Air Pressure at sea level | kPa | MERRA-2, GEOS-5 |
| Specific Humidity | kg/kg | MERRA-2, GEOS-5 |
| Air Temperature | °C | MERRA-2, GEOS-5 |
| Windspeed | m/s | MERRA-2, GEOS-5 |
| Solar Radiation | W/m2  | MERRA-2 |
| Soil Saturation | - | from pywapor.se_root()

#### Temporal SE_Root Data (instantaneous)
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Normalized Difference Vegetation Index (NDVI) | - | MOD13, MYD13, PROBA-V |
| Air Pressure at sea level | kPa | MERRA-2, GEOS-5 |
| Air Pressure at surface level  | kPa | MERRA-2, GEOS-5 |
| Specific Humidity | kg/kg | MERRA-2, GEOS-5 |
| Air Temperature  | °C | MERRA-2, GEOS-5 |
| Windspeed | m/s | MERRA-2, GEOS-5 |
| Total Precipitable Water Vapour  | mm | MERRA-2, GEOS-5 |
| Land Surface Temperature (LST) | K | MOD11, MYD11 |

#### Static Data
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
Landcover | - | WaPOR, GlobCover
Digital Elevation Model | m.a.s.l | SRTM
Air Temperature (yearly amplitude) | K | GLDAS
Latitude | DD | from NDVI
Longitude | DD | from NDVI
Slope | ° | from Elevation
Slope Aspect | ° | from Elevation
Bulk Stomatal Resistance | s/m | from Landcover
Landmask | - | from Landcover
Maximum Light Use Efficiency | gr/MJ | from Landcover
Maximum Obstacle Height | m | from Landcover 

#### 🛰️ Sources
| Source | Temporal Availability | Temporal Resolution |Spatial Resolution | Used For |
| ------ | ------ | ------ | ------ | ------ |
|[MOD13](https://lpdaac.usgs.gov/products/mod13q1v006/) | 2000-02-18 - ongoing | 16-Daily |250m|NDVI|
|[MYD13](https://lpdaac.usgs.gov/products/myd13q1v006/) | 2002-07-04 - ongoing | 16-Daily |250m|NDVI|
|[MCD43](https://lpdaac.usgs.gov/products/mcd43a3v006/)|2000-02-16 - ongoing|Daily|500m|Albedo|
|[MOD11](https://lpdaac.usgs.gov/products/mod11a1v006/) | 2000-02-24 - ongoing | Daily | 1000m | LST |
|[MYD11](https://lpdaac.usgs.gov/products/myd11a1v006/)| 2002-07-04 - ongoing | Daily | 1000m | LST |
|[PROBAV](https://www.vito-eodata.be/collectioncatalogue/srv/eng/catalog.search#/metadata/urn:ogc:def:EOP:VITO:PROBAV_S5-TOC_100M_V001)|2014-03-11 - ongoing|5-Daily|100m|NDVI, Albedo|
| [GEOS5](https://geos5.org) | 2017-12-01 - ongoing | 3-Hourly |0.3125°×0.25° | Meteo |
| [MERRA2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) | 1980-01-01 - ongoing | Hourly | 0.625°×0.5° | Meteo | 
| [CHIRPS](https://www.chc.ucsb.edu/data/chirps) |  1981-01-01 - ongoing | Daily | 0.05° | Precipitation |
| [WAPOR](https://wapor.apps.fao.org/catalog/WAPOR_2/1/L1_LCC_A) | 2009 - 2020 | Yearly |250m | Landcover |
| [GLOBCOVER](http://due.esrin.esa.int/page_globcover.php) | 2009 | Single| 250m | Landcover |
| [SRTM](https://srtm.csi.cgiar.org) | 2009 | Single | 90m | DEM |
