from setuptools import setup, find_packages

setup(
    name = 'pywapor_test',
    version = '2.0.2',
    url = 'https://bitbucket.org/cioapps/wapor-et-look/src/master/',
    author = "FAO",
    author_email = "bert.coerver@fao.org",
    packages = find_packages(include = ['pywapor', 'pywapor.*']),
    install_requires = [
        'gdal>=2.2.0',
        'watertools>=0.0.26',
        'numpy',
        'pandas',
    ]
)