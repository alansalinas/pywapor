from setuptools import setup, find_packages

setup(
    name = 'pywapor',
    version = '2.2.0',
    url = 'https://bitbucket.org/cioapps/wapor-et-look/src/master/',
    author = "FAO",
    author_email = "bert.coerver@fao.org",
    license = "Apache",
    packages = find_packages(include = ['pywapor', 'pywapor.*']),
    include_package_data=True,
    install_requires = [
        'gdal<=3.1.4',
        'numpy',
        'pandas',
        'requests',
        'matplotlib',
        'netcdf4',
        'pyproj',
        'scipy',
        'fiona',
        'pycurl',
        'pyshp',
        'joblib',
        'bs4',
        'paramiko',
        'rasterio',
        'xarray',
        'geojson',
        'vito_download',
        'nest_asyncio',
        'tqdm',
    ]
)