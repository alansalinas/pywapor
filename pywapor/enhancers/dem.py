import numpy as np
import pywapor.general.processing_functions as pf
import xarray as xr

def to_slope(ds, var, out_var = None):
    """Calculate the gradient of `var` (which would usually be `z` or elevation).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing var.
    var : str
        Variable name.
    out_var : str, optional
        New variable name in case `var` should not be overwritten, by default None.

    Returns
    -------
    xr.Dataset
        Enhanced dataset.
    """
    dlat, dlon = pf.calc_dlat_dlon(None, None, None, (ds.lat.values, ds.lon.values))            

    dem = ds[var].values

    x, y = np.gradient(dem, np.nanmean(dlon), np.nanmean(dlat))
    hypotenuse_array = np.hypot(x,y)
    slope = np.arctan(hypotenuse_array) * 180.0 / np.pi

    new_data = xr.DataArray(slope, 
                                    coords = {
                                                "lat": ds.lat, 
                                                "lon": ds.lon})
    
    new_data.attrs = {"sources": ["dem.to_slope()"]}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def to_aspect(ds, var, out_var = None):
    """Calculate the aspect of `var` (which would usually be `z` or elevation).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing var.
    var : str
        Variable name.
    out_var : str, optional
        New variable name in case `var` should not be overwritten, by default None.

    Returns
    -------
    xr.Dataset
        Enhanced dataset.
    """
    dlat, dlon = pf.calc_dlat_dlon(None, None, None, (ds.lat.values, ds.lon.values))            

    dem = ds[var].values

    x, y = np.gradient(dem, np.nanmean(dlon), np.nanmean(dlat))

    aspect = np.arctan2(y/np.nanmean(dlat), 
                        -x/np.nanmean(dlon)) * 180.0 / np.pi + 180.0
    
    new_data = xr.DataArray(aspect, 
                                    coords = {
                                                "lat": ds.lat, 
                                                "lon": ds.lon})
    
    new_data.attrs = {"sources": ["dem.to_aspect()"]}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def to_lat(ds, var, out_var = None):
    """Calculate the latitude in degrees of every pixel.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing var.
    var : str
        Variable name.
    out_var : str, optional
        New variable name in case `var` should not be overwritten, by default None.

    Returns
    -------
    xr.Dataset
        Enhanced dataset.
    """
    lat_deg = np.rot90(np.tile(ds.lat, ds.lon.size).reshape((ds.lon.size,ds.lat.size)), 3)

    new_data = xr.DataArray(lat_deg, 
                                    coords = {
                                                "lat": ds.lat, 
                                                "lon": ds.lon})
    
    new_data.attrs = {"sources": ["dem.to_lat()"]}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data
    
    return ds

def to_lon(ds, var, out_var = None):
    """Calculate the longitude in degrees of every pixel.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing var.
    var : str
        Variable name.
    out_var : str, optional
        New variable name in case `var` should not be overwritten, by default None.

    Returns
    -------
    xr.Dataset
        Enhanced dataset.
    """
    lon_deg = np.tile(ds.lon, ds.lat.size).reshape((ds.lat.size,ds.lon.size))

    new_data = xr.DataArray(lon_deg, 
                                    coords = {
                                                "lat": ds.lat, 
                                                "lon": ds.lon})
    
    new_data.attrs = {"sources": ["dem.to_lon()"]}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds