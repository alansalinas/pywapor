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