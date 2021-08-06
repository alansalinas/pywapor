from setuptools import setup, find_packages

setup(
    name = 'pywapor',
    version = '2.1.1',
    url = 'https://bitbucket.org/cioapps/wapor-et-look/src/master/',
    author = "FAO",
    author_email = "bert.coerver@fao.org",
    license = "Apache v2.0",
    packages = find_packages(include = ['pywapor', 'pywapor.*']),
    install_requires = [
        'gdal<=3.1.4',
        'watertools<=0.0.28',
        'numpy',
        'pandas',
        'requests',
        'matplotlib',
        # These should be requirements for watertools, not pywapor:
        'netcdf4',
        'pyproj',
        'scipy',
        'fiona',
        'pycurl',
        'pyshp',
        'joblib',
        'bs4',
        'paramiko',
    ]
)