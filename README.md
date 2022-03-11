## pyWaPOR

![downloads](https://img.shields.io/pypi/dw/pywapor) [![version](https://img.shields.io/pypi/v/pywapor)](https://pypi.org/project/pywapor/) ![python](https://img.shields.io/pypi/pyversions/pywapor) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/1_introduction.ipynb)  

This repository contains a Python implementation of the algorithm used to generate the [WaPOR](http://www.fao.org/in-action/remote-sensing-for-water-productivity/en/) [datasets](https://wapor.apps.fao.org/home/WAPOR_2/1). It can be used to calculate evaporation, transpiration and biomass production maps.

### Installation

Its recommended to install in a clean [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and use [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install GDAL before installing pywapor. Only Python 3.7 and 3.8 are supported.

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

### Usage

To run the model for one dekad (from 2021-07-01 to 2021-07-11 in this case) for the Fayoum irrigation scheme in Egypt (but feel free to change the [boundingbox](http://bboxfinder.com) defined by `latlim` and `lonlim`) using mainly MODIS data, run the following code. 

```python
import pywapor

# User inputs.
startdate = "2021-07-01"
enddate = "2021-07-11"
latlim = [28.9, 29.7]
lonlim = [30.2, 31.2]
project_folder = r"/my_first_ETLook_run/"

# Download and prepare input data.
ds_in, fh_in = pywapor.pre_et_look.main(project_folder, startdate, enddate, latlim, lonlim)

# Run the model.
ds_out = pywapor.et_look.main(ds_in)
```

Check out the documentation and the notebooks below to learn more!

### Documentation
Go [here](https://www.fao.org/aquastat/py-wapor/) for the full pyWaPOR documentation.

#### Notebooks

<table class = "docutils align-default">
   <thead>
      <tr class="row-odd" style="text-align:center">
         <th class="head"></th>
         <th class="head">Name</th>
         <th class="head">Duration<sup>1</sup></th>
         <th class="head" width = "150">Colab</th>
      </tr>
   </thead>
   <tbody>
      <tr class="row-odd">
         <td>1.</td>
         <td>Introduction</td>
         <td>10 + 100</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/1_introduction.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
      <tr class="row-even">
         <td>2.</td>
         <td>Composites</td>
         <td>10 + 10</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/2_composites.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
      <tr class="row-odd">
         <td>3.</td>
         <td>Levels</td>
         <td>10 + 20</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/3_levels.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
      <tr class="row-even">
         <td>4.</td>
         <td>Sideloading</td>
         <td>10 + 5</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/4_sideloading.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
      <tr class="row-odd">
         <td>5.</td>
         <td>Enhancers</td>
         <td>10 + 5</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/5_enhancers.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
      <tr class="row-even">
         <td>6.</td>
         <td>pyWaPOR vs. WaPOR</td>
         <td>10 + 10</td>
         <td style="text-align:center"><a href="https://colab.research.google.com/github/bertcoerver/pywapor_notebooks/blob/main/6_wapor_vs_pywapor.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="colab"/></a></td>
      </tr>
   </tbody>
</table>

<sup>1</sup> Estimation of the time required in minutes, as in "active" + "download time", assuming the notebooks have access to previously downloaded data.

#### WaPOR v2

<ul>
<li><a href="https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_ET_v2_data_manual_finaldraft_v2.2.pdf">WaPOR-ETLook Data Manual (v2)</a></li>

<li><a href="https://bitbucket.org/cioapps/wapor-et-look/downloads/FRAME_NPP_v2_data_manual_finaldraft_v2.2.pdf">WaPOR-Biomass Data Manual (v2)</a></li>
</ul>

#### WaPOR v1

<ul>
<li><a href="https://bitbucket.org/cioapps/wapor-et-look/raw/9ec88e56769f49722c2d1165bb34547f5842b811/Docs/WaPOR_ET_data_manual_finaldraft-v1.2-for-distribution.pdf">WaPOR-ETLook Data Manual (v1)</a></li>
</ul>

### Acknowledgments
The methodology for WaPOR was developed by the FRAME1 consortium, consisting of eLEAF (lead), VITO, ITC, University of Twente and Waterwatch foundation, commissioned by and in partnership with the Land and Water Division of FAO. The method for calculating evapotranspiration is based on the ETLook model developed by eLEAF in 2010. The method for calculating total biomass production is based on the C-Fix model. 

The code in the pywapor.et_look_v2 module of this repository, containing all core physical functions used by ETLook, was written by Henk Pelgrum (eLEAF) and Rutger Kassies (eLEAF). The remaining modules have been developed by Bert Coerver (FAO), Tim Hessels (WaterSat), and, in the framework of the ESA-funded ET4FAO project, Radoslaw Guzinski (DHI-GRAS), Hector Nieto (Complutig) and Laust Faerch (DHI-GRAS).

### Contact
For questions, requests or issues with this repository, please contact Bert Coerver at [bert.coerver@fao.org](mailto:bert.coerver@fao.org) or the WaPOR team at [wapor@fao.org](mailto:wapor@fao.org).

### Release Notes

#### 2.4.1 (2022-03-11)
<br>
<ul>
    <li> NetCDF files are now compressed when saved to disk.</li>
    <li> Calculation of Total Biomass Production is now turned on by default.</li>
    <li> It is no longer required to provide <b>all</b> input variables to et_look,
    the model will calculate as many variables as possible with the given data. For example,
    if you are only interested in acquiring interception rates, it would now suffice to only prepare ndvi and p (precipitation) data with pre_et_look.</li>
    <li> et_look now automatically generates an interactive network graph visualising the executed computation steps.</li>
</ul>

#### 2.4.0 (2022-02-03)
<br>
<ul>
    <li> Easily apply your own functions to data, i.e. use your own custom filters, gap-fillers etc.</li>
    <li> Side-load your own data, i.e. easily incorporate you own datasets.</li>
    <li> Added functions to process Landsat Level-2 to “ndvi”, “lst” and “r0”.</li>
    <li> Data is now stored and processed as netCDF (using xarray and dask).</li>
    <li> Calculations in et_look() and se_root() are now done in chunks, instead of using a for-loop.</li>
    <li> Some previously constant parameters now have spatial variability.</li>
    <li> Improved logging.</li>
    <li> Download functions now show progress and download-speed.</li>
    <li> MODIS datasets switched from v6.0 (decommissioned soon) to v6.1.</li>
    <li> The lapse-rate correction to temperature data is now more accurate and faster.</li>
    <li> VITO and WAPOR passwords are now checked when entered.</li>
    <li> Other bug-fixes and performance improvements.</li>
</ul>

#### 2.3.0 (2021-11-19)
<br>
<ul> 
    <li>Automatically create input composites before running ETLook.</li>
    <li>Choose composite lengths in number of days or dekads.</li>
    <li>Option to choose which products to use per variable.</li>
    <li>Calculate soil saturation separate from ETLook.</li>
    <li>PROBA-V support for NDVI and Albedo inputs.</li>
    <li>Define diagnostics pixels, for which extra outputs are created (e.g. charts, maps etc.).</li>
    <li>Bug-fixes and performance improvements.</li>
</ul>

### License
[APACHE](https://bitbucket.org/cioapps/wapor-et-look/src/dev/LICENSE)