import numpy as np
from scipy.signal import convolve2d
import numpy as np
import xarray as xr
from pywapor.general.logger import log

def kelvin_to_celsius(ds, var, out_var = None):
    """Convert units of a variable from Kelvin to degrees Celsius.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the data.
    var : str
        Name of the variable inside ds to apply this function to.
    out_var : str, optional
        Instead of overwriting `var`, store the new data in `out_var`, 
        by default None.

    Returns
    -------
    xarray.Dataset
        Adjusted data.
    """

    new_data = ds[var] - 273.15
    new_data.attrs = {"unit": "°C"}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def celsius_to_kelvin(ds, var, out_var = None):
    """Convert units of a variable from Kelvin to degrees Celsius.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the data.
    var : str
        Name of the variable inside ds to apply this function to.
    out_var : str, optional
        Instead of overwriting `var`, store the new data in `out_var`, 
        by default None.

    Returns
    -------
    xarray.Dataset
        Adjusted data.
    """

    new_data = ds[var] + 273.15
    new_data.attrs = {"unit": "K"}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def lapse_rate(ds, var, out_var = None, lapse_var = "z", 
                radius = 0.25, lapse_rate = -0.006):
    """Apply a lapse-rate to `var` inside `ds` based on the difference between
    the actual value of `lapse_var` and its local (defined by `radius`) mean.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with the data.
    var : str
        Name of the variable inside ds to apply this function to, i.e. y.
    out_var : str, optional
        Instead of overwriting `var`, store the new data in `out_var`, 
        by default None.
    lapse_var : str, optional
        The variable used to calculate dx.
    radius : float, optional
        The local altitude is compared to the local mean altitude, radius
        (in degrees) defines over which area to calculate the mean, by default 0.25.
    lapse_rate : float, optional
        dy/dx, i.e. how much should `var` be adjusted per difference in 
        `lapse_var`, by default -0.006.

    Returns
    -------
    xarray.Dataset
        Adjusted data.
    """
    
    if "t_diff" not in list(ds.keys()):

        log.info(f"--> Calculating local means (r = {radius:.2f}°) of `{lapse_var}`.")
        
        pixel_size = (np.median(np.diff(ds.lon)), 
                        np.median(np.diff(ds.lat)))
        dem = ds[lapse_var].values
        
        dem_coarse = local_mean(dem, pixel_size, radius)
        t_diff = (dem - dem_coarse) * lapse_rate
        ds["t_diff"] = xr.DataArray(t_diff, 
                                    coords = {
                                                "lat": ds.lat, 
                                                "lon": ds.lon})

    new_data = ds[var] + ds["t_diff"]

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else:
        ds[var] = new_data

    return ds

def create_circle(pixel_size_x, pixel_size_y, radius, boolean = False):
    """Create a 2d array with pixel values representing the distance to 
    the center of the array.

    Parameters
    ----------
    pixel_size_x : float
        Size of a pixel in the x-direction. Unit can be anything, as long as it
        is the same as for `pixel_size_y` and `radius`.
    pixel_size_y : float
        Size of a pixel in the y-direction. Unit can be anything, as long as it
        is the same as for `pixel_size_x` and `radius`.
    radius : float
        Size of the circle. Unit can be anything, as long as it
        is the same as for `pixel_size_x` and `pixel_size_y`.
    boolean : bool, optional
        Return distances to center (False) or whether or not the distance to 
        the center is smaller than `radius` (True), by default False.

    Returns
    -------
    ndarray
        Values representing the distance to the center of the array.
    """
    min_pixel_size = np.min([np.abs(pixel_size_y), np.abs(pixel_size_x)])

    length = int(np.ceil(radius / min_pixel_size)) + 1

    Y, X = np.ogrid[:length, :length]  * np.array([pixel_size_x, pixel_size_y])

    center = (Y[-1, -1], X[-1,-1])

    dist = np.sqrt((X - center[1])**2 + (Y-center[0])**2)

    half_circle = np.concatenate([dist, np.flipud(dist)[1:,:]])
    circle = np.concatenate([half_circle, np.fliplr(half_circle)[:,1:]], axis =1)
    
    if boolean:
        circle = (circle <= radius).astype(int)
    
    return circle

def local_mean(array, pixel_size, radius):
    """Calculate local means.

    Parameters
    ----------
    array : ndarray
        Two dimensional array for which local means will be calculated.
    pixel_size : tuple
        Size of the pixels in `array` in x and y direction respectively.
    radius : float
        Size of the circle around each pixel from which to 
        calculate the local mean.

    Returns
    -------
    ndarray :
        Array with the local mean of each pixel.
    """

    pixel_size_x, pixel_size_y = pixel_size
    circle_weights = create_circle(pixel_size_x, pixel_size_y, 
                                    radius, boolean=True)
    
    # if isinstance(name, type(None)):
    #     log.info("--> Calculating local means (step 1 of 2).")
    # else:
    #     log.info(f"--> Calculating local means of `{name}` (step 1/2).")
    divisor = convolve2d(np.ones_like(array), circle_weights, mode = "same")
    
    # if isinstance(name, type(None)):
    #     log.info("--> Calculating local means (step 2 of 2).")
    # else:
    #     log.info(f"--> Calculating local means of `{name}` (step 2/2).") 
    array_coarse = convolve2d(array, circle_weights, mode = "same") / divisor
    return array_coarse
