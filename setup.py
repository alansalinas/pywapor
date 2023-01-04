from setuptools import setup, find_packages

setup(
    name = 'pywapor',
    version = '3.2.0',
    url = 'https://www.fao.org/aquastat/py-wapor/',
    author = "FAO",
    author_email = "bert.coerver@fao.org",
    license = "Apache",
    packages = find_packages(include = ['pywapor', 'pywapor.*']),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires = [
        'gdal',
        'xarray>=0.20',
        'numpy',
        'pydap',
        'pandas',
        'requests',
        'matplotlib',
# NOTE otherwise opendap gives problem in colab, in conda env netcdf=1.6.0 
# works fine -> https://github.com/Unidata/netcdf4-python/issues/1179
# NOTE also, cant install netcdf4 with conda, because the conda-forge distribution cant open
# the PROBAV HDF files. See issues https://github.com/Unidata/netcdf4-python/issues/1182 and
# https://github.com/conda-forge/netcdf4-feedstock/issues/136
        'netcdf4<1.6.0', 
        'pyproj',
        'scipy',
        'pycurl',
        'pyshp',
        'joblib',
        'bs4',
        'rasterio',
        'bottleneck>=1.3.1',
        'geojson',
        'tqdm',
        'dask',
        'rioxarray',
        'python_log_indenter',
        'cryptography',
        'pyvis',
        'cachetools',
        'cdsapi',
        'sentinelsat',
        'shapely',
        'lxml',
        'geopy',
        'scikit-learn',
        'numba',
        'xmltodict',
# NOTE Another fix for Colab... https://github.com/googlecolab/colabtools/issues/3134
        'importlib-metadata==4.13.0',
    ],
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)