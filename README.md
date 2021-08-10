# pyWAPOR
![downloads](https://img.shields.io/pypi/dw/pywapor)
[![version](https://img.shields.io/pypi/v/pywapor)](https://pypi.org/project/pywapor/)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1YEsCN6GnMGvOzXT4YaJ_jeu58mGIhhMq?usp=sharing)


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

To run the model for one day (2019-07-07 in this case) for the Fayoum irrigation scheme in Egypt (but feel free to change the [boundingbox](http://bboxfinder.com) defined by `latlim` and `lonlim`) using MODIS data, use the following code. 

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

See the `examples` folder for more examples or check out the [Colab Notebook](https://colab.research.google.com/drive/1YEsCN6GnMGvOzXT4YaJ_jeu58mGIhhMq?usp=sharing).

## Data Sources

To run the ETLook model, two types of spatial data are required, temporal and static data. **Each of these variables can be collected from whichever source you wish to use**, as long as you make sure the units are correct, the data is stored as a GeoTIFF (1 band per file, 1 file for each variable and date), the files all have the same no-data-value and they all have the same projection and resolution. The `tests/test_data/input/` folder contains an example input dataset.

**For your convenience, the pyWAPOR package has a function that can collect all this data from selected sources** and make sure the data is stored in the correct format and folder structure.

#### Temporal Data
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Albedo | - | MODIS, *️⃣ |
| Land Surface Temperature | K | MODIS, *️⃣ |
| Normalized Difference Vegetation Index | - | MODIS, *️⃣ |
| Air Pressure at sea level (daily average) | kPa | MERRA-2, GEOS-5 |
| Air Pressure at sea level (instanteneous) | kPa | MERRA-2, GEOS-5 |
| Air Pressure at surface level (instanteneous) | kPa | MERRA-2, GEOS-5 |
| Precipitation | mm/day | CHIRPS |
| Specific Humidity (daily average) | kg/kg | MERRA-2, GEOS-5 |
| Specific Humidity (instanteneous) | kg/kg | MERRA-2, GEOS-5 |
| Air Temperature (daily average) | °C | MERRA-2, GEOS-5 |
| Air Temperature (instanteneous) | °C | MERRA-2, GEOS-5 |
| Air Temperature (daily maximum) | °C | MERRA-2, GEOS-5 |
| Air Temperature (daily minimum) | °C | MERRA-2, GEOS-5 |
| Transmissivity | - | MERRA-2 |
| Windspeed (daily average) | m/s | MERRA-2, GEOS-5 |
| Windspeed (instanteneous) | m/s | MERRA-2, GEOS-5 |
| Total Precipitable Water Vapout | mm | MERRA-2, GEOS-5 |
| Instantaneous Data Time | hour | n/a

*️⃣ PROBA-V and LandSat support in development.

#### Static Data
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
Landcover | - | WaPOR, GlobCover
Digital Elevation Model | m.a.s.l | SRTM
Air Temperature (yearly amplitude) | K | MERRA-2
Latitude | DD | from Albedo
Longitude | DD | from Albedo
Slope | ° | from Elevation
Slope Aspect | ° | from Elevation
Bulk Stomatal Resistance | s/m | from Landcover
Landmask | - | from Landcover
Maximum Light Use Efficiency | gr/MJ | from Landcover
Maximum Obstacle Height | m | from Landcover

#### 🛰️ Sources
| Source | Temporal Availability | Spatial Resolution | Used For |
| ------ | ------ | ------ | ------ |
| [MODIS](https://modis.gsfc.nasa.gov) | 2000-02-18 - ongoing | ~250m, ~500m, ~1000m | NDVI, Albedo, LST |
| [GEOS-5](https://geos5.org) | 2017-12-01 - ongoing | 0.3125°×0.25° | Meteo |
| [MERRA-2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) | 1980-01-01 - ongoing | 0.625°×0.5° | Meteo | 
| [CHIRPS](https://www.chc.ucsb.edu/data/chirps) |  1981-01-01 - ongoing | 0.05°×0.05° | Precipitation |
| [WaPOR](https://wapor.apps.fao.org/catalog/WAPOR_2/1/L1_LCC_A) | 2009 - 2020 | ~250m | Landcover |
| [GlobCover](http://due.esrin.esa.int/page_globcover.php) | 2009 | ~250m | Landcover |
| [SRTM](https://srtm.csi.cgiar.org) | 2009 | ~90m | DEM |

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
