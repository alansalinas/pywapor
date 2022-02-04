"""Functions and dictionaries to convert Landsat Quality Assesment Bitmasks
into numpy boolean arrays.
"""

import matplotlib.pyplot as plt
import numpy as np

def get_pixel_qa_bits(collection, ls_number, level):
    """Returns a dictionary to bridge between labels and bits in
    Landsat pixel_qa bands.

    Parameters
    ----------
    collection : int
        For which collection to get the values.
    ls_number : int
        For which spacecraft to get the values.
    level : int
        For which level to get the values.

    Returns
    -------
    dict
        Keys are labels (different per Landsat mission/product) and values are
        definitions of the bits.

    Examples
    --------
    # Count from right to left.
    # <-etc...9876543210
    #

    Returns True if 2nd bit is set:
    [(0b0000000000000100, False)]

    Returns True if 9th bit is set AND 8th bit is unset.
    [(0b0000001000000000, False)
     (0b0000000100000000, True)]

    Returns True if 11th bit is unset AND 10th bit is unset.
    [(0b0000100000000000, True)
     (0b0000010000000000, True)]
    """


    c2l8lvl1 = {
        'fill':             [(0b0000000000000001, False)],
        'dilated_cloud':    [(0b0000000000000010, False)],
        'cirrus':           [(0b0000000000000100, False)],
        'cloud':            [(0b0000000000001000, False)],
        'cloud_shadow':     [(0b0000000000010000, False)],
        'snow':             [(0b0000000000100000, False)], 
        'clear':            [(0b0000000001000000, False)],
        'water':            [(0b0000000010000000, False)],
        'cloud_nd':         [(0b0000001000000000, True),
                             (0b0000000100000000, True)],
        'cloud_low':        [(0b0000001000000000, True),
                             (0b0000000100000000, False)],
        'cloud_med':        [(0b0000001000000000, False),
                             (0b0000000100000000, True)],
        'cloud_high':       [(0b0000001000000000, False),
                             (0b0000000100000000, False)],
        'cloud_shadow_nd':  [(0b0000100000000000, True),    
                             (0b0000010000000000, True)],   
        'cloud_shadow_low': [(0b0000100000000000, True),    
                             (0b0000010000000000, False)],  
        'cloud_shadow_med': [(0b0000100000000000, False),   
                             (0b0000010000000000, True)],   
        'cloud_shadow_high':[(0b0000100000000000, False),   
                             (0b0000010000000000, False)],  
        'snow_ice_nd':      [(0b0010000000000000, True),    
                             (0b0001000000000000, True)],   
        'snow_ice_low':     [(0b0010000000000000, True),    
                             (0b0001000000000000, False)],  
        'snow_ice_med':     [(0b0010000000000000, False),   
                             (0b0001000000000000, True)],   
        'snow_ice_high':    [(0b0010000000000000, False),   
                             (0b0001000000000000, False)],  
        'cirrus_nd':        [(0b1000000000000000, True),    
                             (0b0100000000000000, True)],   
        'cirrus_low':       [(0b1000000000000000, True),    
                             (0b0100000000000000, False)],  
        'cirrus_med':       [(0b1000000000000000, False),   
                             (0b0100000000000000, True)],   
        'cirrus_high':      [(0b1000000000000000, False),   
                             (0b0100000000000000, False)],  
    }

    c2l7lvl1 = {
        'fill':             [(0b0000000000000001, False)],
        'dilated_cloud':    [(0b0000000000000010, False)],
        'cloud':            [(0b0000000000001000, False)],
        'cloud_shadow':     [(0b0000000000010000, False)],
        'snow':             [(0b0000000000100000, False)], 
        'clear':            [(0b0000000001000000, False)],
        'water':            [(0b0000000010000000, False)],
        'cloud_nd':         [(0b0000001000000000, True),
                             (0b0000000100000000, True)],
        'cloud_low':        [(0b0000001000000000, True),
                             (0b0000000100000000, False)],
        'cloud_med':        [(0b0000001000000000, False),
                             (0b0000000100000000, True)],
        'cloud_high':       [(0b0000001000000000, False),
                             (0b0000000100000000, False)],
        'cloud_shadow_nd':  [(0b0000100000000000, True),    
                             (0b0000010000000000, True)],   
        'cloud_shadow_low': [(0b0000100000000000, True),    
                             (0b0000010000000000, False)],  
        'cloud_shadow_med': [(0b0000100000000000, False),   
                             (0b0000010000000000, True)],   
        'cloud_shadow_high':[(0b0000100000000000, False),   
                             (0b0000010000000000, False)],  
        'snow_ice_nd':      [(0b0010000000000000, True),    
                             (0b0001000000000000, True)],   
        'snow_ice_low':     [(0b0010000000000000, True),    
                             (0b0001000000000000, False)],  
        'snow_ice_med':     [(0b0010000000000000, False),   
                             (0b0001000000000000, True)],   
        'snow_ice_high':    [(0b0010000000000000, False),   
                             (0b0001000000000000, False)],   
    }

    c1l8lvl1 = {
        'designated_fill':  [(0b0000000000000001, False)],  
        'terrain_occlusion':[(0b0000000000000010, False)],  
        'rad_sat_0':        [(0b0000000000001000, True),    
                             (0b0000000000000100, True)],   
        'rad_sat_12':       [(0b0000000000001000, True),    
                             (0b0000000000000100, False)],  
        'rad_sat_34':       [(0b0000000000001000, False),   
                             (0b0000000000000100, True)],   
        'rad_sat_5':        [(0b0000000000001000, False),   
                             (0b0000000000000100, False)],  
        'cloud':            [(0b0000000000010000, False)],  
        'cloud_nd':         [(0b0000000001000000, True),    
                             (0b0000000000100000, True)],   
        'cloud_low':        [(0b0000000001000000, True),    
                             (0b0000000000100000, False)],  
        'cloud_med':        [(0b0000000001000000, False),   
                             (0b0000000000100000, True)],   
        'cloud_high':       [(0b0000000001000000, False),   
                             (0b0000000000100000, False)],  
        'cloud_shadow_nd':  [(0b0000000100000000, True),    
                             (0b0000000010000000, True)],   
        'cloud_shadow_low': [(0b0000000100000000, True),    
                             (0b0000000010000000, False)],  
        'cloud_shadow_med': [(0b0000000100000000, False),   
                             (0b0000000010000000, True)],   
        'cloud_shadow_high':[(0b0000000100000000, False),   
                             (0b0000000010000000, False)],  
        'snow_ice_nd':      [(0b0000010000000000, True),    
                             (0b0000001000000000, True)],   
        'snow_ice_low':     [(0b0000010000000000, True),    
                             (0b0000001000000000, False)],  
        'snow_ice_med':     [(0b0000010000000000, False),   
                             (0b0000001000000000, True)],   
        'snow_ice_high':    [(0b0000010000000000, False),   
                             (0b0000001000000000, False)],  
        'cirrus_nd':        [(0b0001000000000000, True),    
                             (0b0000100000000000, True)],   
        'cirrus_low':       [(0b0001000000000000, True),    
                             (0b0000100000000000, False)],  
        'cirrus_med':       [(0b0001000000000000, False),   
                             (0b0000100000000000, True)],   
        'cirrus_high':      [(0b0001000000000000, False),   
                             (0b0000100000000000, False)],
    }

    c1l754lvl1 = {
        'designated_fill':  [(0b0000000000000001, False)],  
        'designated_pixel': [(0b0000000000000010, False)],  
        'rad_sat_0':        [(0b0000000000001000, True),    
                             (0b0000000000000100, True)],   
        'rad_sat_12':       [(0b0000000000001000, True),    
                             (0b0000000000000100, False)],  
        'rad_sat_34':       [(0b0000000000001000, False),   
                             (0b0000000000000100, True)],   
        'rad_sat_5':        [(0b0000000000001000, False),   
                             (0b0000000000000100, False)],  
        'cloud':            [(0b0000000000010000, False)],  
        'cloud_nd':         [(0b0000000001000000, True),    
                             (0b0000000000100000, True)],   
        'cloud_low':        [(0b0000000001000000, True),    
                             (0b0000000000100000, False)],  
        'cloud_med':        [(0b0000000001000000, False),   
                             (0b0000000000100000, True)],   
        'cloud_high':       [(0b0000000001000000, False),   
                             (0b0000000000100000, False)],  
        'cloud_shadow_nd':  [(0b0000000100000000, True),    
                             (0b0000000010000000, True)],   
        'cloud_shadow_low': [(0b0000000100000000, True),    
                             (0b0000000010000000, False)],  
        'cloud_shadow_med': [(0b0000000100000000, False),   
                             (0b0000000010000000, True)],   
        'cloud_shadow_high':[(0b0000000100000000, False),   
                             (0b0000000010000000, False)],  
        'snow_ice_nd':      [(0b0000010000000000, True),    
                             (0b0000001000000000, True)],   
        'snow_ice_low':     [(0b0000010000000000, True),    
                             (0b0000001000000000, False)],  
        'snow_ice_med':     [(0b0000010000000000, False),   
                             (0b0000001000000000, True)],   
        'snow_ice_high':    [(0b0000010000000000, False),   
                             (0b0000001000000000, False)],  
    }

    c2l8lvl2 = {
        'fill':             [(0b0000000000000001, False)],  
        'dilated_cloud':    [(0b0000000000000010, False)],  
        'cirrus':           [(0b0000000000000100, False)],   
        'cloud':            [(0b0000000000001000, False)],  
        'cloud_shadow':     [(0b0000000000010000, False)],   
        'snow':             [(0b0000000000100000, False)],  
        'clear':            [(0b0000000001000000, False)],  
        'water':            [(0b0000000010000000, False)],
  
        'cloud_nd':         [(0b0000000100000000, False),    
                             (0b0000001000000000, False)],  
        'cloud_low':        [(0b0000000100000000, True),   
                             (0b0000001000000000, False)],   
        'cloud_med':        [(0b0000000100000000, False),   
                             (0b0000001000000000, True)],  
        'cloud_high':       [(0b0000000100000000, True),    
                             (0b0000001000000000, True)],

        'cloud_shadow_nd':  [(0b0000010000000000, False),    
                             (0b0000100000000000, False)],  
        'cloud_shadow_low': [(0b0000010000000000, True),   
                             (0b0000100000000000, False)],   
        'cloud_shadow_med': [(0b0000010000000000, False),   
                             (0b0000100000000000, True)],  
        'cloud_shadow_high':[(0b0000010000000000, True),    
                             (0b0000100000000000, True)],   

        'snow_ice_nd':      [(0b0001000000000000, False),    
                             (0b0010000000000000, False)],  
        'snow_ice_low':     [(0b0001000000000000, True),   
                             (0b0010000000000000, False)],   
        'snow_ice_med':     [(0b0001000000000000, False),   
                             (0b0010000000000000, True)],
        'snow_ice_high':    [(0b0001000000000000, True),    
                             (0b0010000000000000, True)],

        'cirrus_nd':        [(0b0100000000000000, False),   
                             (0b1000000000000000, False)],   
        'cirrus_low':       [(0b0100000000000000, True),   
                             (0b1000000000000000, False)], 
        'cirrus_med':       [(0b0100000000000000, False),   
                             (0b1000000000000000, True)],   
        'cirrus_high':      [(0b0100000000000000, True),   
                             (0b1000000000000000, True)],  
    }

    c2l7lvl2 = {
        'fill':             [(0b0000000000000001, False)],  
        'dilated_cloud':    [(0b0000000000000010, False)],    
        'cloud':            [(0b0000000000001000, False)],  
        'cloud_shadow':     [(0b0000000000010000, False)],   
        'snow':             [(0b0000000000100000, False)],  
        'clear':            [(0b0000000001000000, False)],  
        'water':            [(0b0000000010000000, False)],
  
        'cloud_nd':         [(0b0000000100000000, False),    
                             (0b0000001000000000, False)],  
        'cloud_low':        [(0b0000000100000000, True),   
                             (0b0000001000000000, False)],   
        'cloud_med':        [(0b0000000100000000, False),   
                             (0b0000001000000000, True)],  
        'cloud_high':       [(0b0000000100000000, True),    
                             (0b0000001000000000, True)],

        'cloud_shadow_nd':  [(0b0000010000000000, False),    
                             (0b0000100000000000, False)],  
        'cloud_shadow_low': [(0b0000010000000000, True),   
                             (0b0000100000000000, False)],   
        'cloud_shadow_med': [(0b0000010000000000, False),   
                             (0b0000100000000000, True)],  
        'cloud_shadow_high':[(0b0000010000000000, True),    
                             (0b0000100000000000, True)],   

        'snow_ice_nd':      [(0b0001000000000000, False),    
                             (0b0010000000000000, False)],  
        'snow_ice_low':     [(0b0001000000000000, True),   
                             (0b0010000000000000, False)],   
        'snow_ice_med':     [(0b0001000000000000, False),   
                             (0b0010000000000000, True)],
        'snow_ice_high':    [(0b0001000000000000, True),    
                             (0b0010000000000000, True)],

    }

    all_flags =dict()
    all_flags[1] = dict()
    all_flags[2] = dict()
    all_flags[1][8] = dict()
    all_flags[2][8] = dict()
    all_flags[2][7] = dict()
    all_flags[1][7] = dict()
    all_flags[1][5] = dict()
    all_flags[1][4] = dict()
    all_flags[1][8][1] = c1l8lvl1
    all_flags[2][8][1] = c2l8lvl1
    all_flags[2][7][1] = c2l7lvl1
    all_flags[1][7][1] = c1l754lvl1
    all_flags[1][5][1] = c1l754lvl1
    all_flags[1][4][1] = c1l754lvl1
    all_flags[2][8][2] = c2l8lvl2
    all_flags[2][7][2] = c2l7lvl2

    return all_flags[collection][ls_number][level]

def get_radsat_qa_bits(collection, ls_number, level):
    """Returns a dictionary to bridge between labels and bits in
    Landsat radsat_qa bands.

    Parameters
    ----------
    collection : int
        For which collection to get the values.
    ls_number : int
        For which spacecraft to get the values.
    level : int
        For which level to get the values.

    Returns
    -------
    dict
        Keys are labels (different per Landsat mission/product) and values are
        definitions of the bits.

    Examples
    --------
    # Count from right to left.
    # <-etc...9876543210
    #

    Returns True if 2nd bit is set:
    [(0b0000000000000100, False)]

    Returns True if 9th bit is set AND 8th bit is unset.
    [(0b0000001000000000, False)
     (0b0000000100000000, True)]

    Returns True if 11th bit is unset AND 10th bit is unset.
    [(0b0000100000000000, True)
     (0b0000010000000000, True)]
    """

    c2l7lvl2 = {
        'saturated_band1'  : [(0b0000000000000001, False)],
        'saturated_band2'  : [(0b0000000000000010, False)],
        'saturated_band3'  : [(0b0000000000000100, False)],
        'saturated_band4'  : [(0b0000000000001000, False)],
        'saturated_band5'  : [(0b0000000000010000, False)],
        'saturated_band6L' : [(0b0000000000100000, False)], 
        'saturated_band7'  : [(0b0000000001000000, False)],
        'saturated_band6H' : [(0b0000000100000000, False)],
        'dropped_pixel'    : [(0b0000001000000000, False)],
    }

    c2l8lvl2 = {
        'saturated_band1'  : [(0b0000000000000001, False)],
        'saturated_band2'  : [(0b0000000000000010, False)],
        'saturated_band3'  : [(0b0000000000000100, False)],
        'saturated_band4'  : [(0b0000000000001000, False)],
        'saturated_band5'  : [(0b0000000000010000, False)],
        'saturated_band6'  : [(0b0000000000100000, False)], 
        'saturated_band7'  : [(0b0000000001000000, False)],
        'saturated_band9'  : [(0b0000000100000000, False)],
        'terrain_occlusion': [(0b0000100000000000, False)],
    }

    all_flags =dict()
    all_flags[2] = dict()
    all_flags[2][8] = dict()
    all_flags[2][7] = dict()
    all_flags[2][7][2] = c2l7lvl2
    all_flags[2][8][2] = c2l8lvl2

    return all_flags[collection][ls_number][level]

def open_array(path, slc = [slice(None), slice(None)]):
    """Open an image as an numpy array.

    Parameters
    ----------
    path : str
        Path to image file.
    slc : list, optional
        Two slice objects to return only a part of the image, by 
        default [slice(None), slice(None)]

    Returns
    -------
    np.ndarray
        The pixel values of the image.
    """
    array = plt.imread(path)[slc[0],slc[1],...]
    return array

def get_mask(qa_array, flags, flag_bits):
    """Given a Landsat bitmask ('qa_array') a list of 'flags' to identify and
    a dictionary to convert between the bits and labels, returns a 
    boolean array.

    Parameters
    ----------
    qa_array : np.ndarray
        The landsat bitmask
    flags : [type]
        Which classes to look for, e.g. ["cloud", "cloud_shadow"].
    flag_bits : dict
        Keys are labels (defining the valid items in 'flags'), values are bits.
        See 'ls_bitmasks.get_pixel_qa_bits()'.

    Returns
    -------
    np.ndarray
        Boolean array with True for pixels belonging to the classes 
        specified in 'flags'.
    """

    final_mask = np.zeros_like(qa_array)

    for flag in flags:
        all_checks = list()
        for byte, inverse in flag_bits[flag]:
            if inverse:
                all_checks.append(np.bitwise_and(~qa_array, byte) > 0 )
            else:
                all_checks.append(np.bitwise_and(qa_array, byte) > 0 )
        flag_mask = np.all(all_checks, axis = 0)
        final_mask = final_mask | flag_mask
    
    return final_mask > 0

def create_masked_rgb(red, green, blue, mask = None):
    """Create a 3-dimensional array from three different bands, optionally
    specifying (transparent) colors to use for masked pixels.

    Parameters
    ----------
    red : np.ndarray
        Pixel values to plot in the red channel.
    green : np.ndarray
        Pixel values to plot in the red channel.
    blue : np.ndarray
        Pixel values to plot in the red channel.
    mask : dictionary, optional
        Should contain; a key `array` with a boolean np.ndarray with the same
        shape as `red`, `blue` and `green`; `alpha` with a float value and; 
        `color` as a tuple with rgb values. By default None

    Returns
    -------
    np.ndarray
        3-dimensional array with pixel values to be plotted in red, green 
        and blue channels.

    Example
    -------
    >>> my_mask = np.array([[True, False],[False, True]])
    >>> mask_plot_info = {"array": my_mask,
                          "color": (1,0,0),
                          "alpha": 0.4}
    >>> red = np.array([[200, 500],[240, 34]])
    >>> green = np.array([[120, 500],[240, 34]])
    >>> blue = np.array([[7856, 5650],[2040, 9486]])
    >>> rgb = create_masked_rgb(red, green, blue, mask = mask_plot_info)
    >>> plt.imshow(rgb)
    """

    R = (red - np.min(red))/np.ptp(red)
    G = (green - np.min(green))/np.ptp(green)
    B = (blue - np.min(blue))/np.ptp(blue)

    if not isinstance(mask, type(None)):

        R = np.where(mask["array"] == True, np.clip(R*(1-mask["alpha"])+mask["color"][0]*mask["alpha"], 0, 1), R)
        G = np.where(mask["array"] == True, np.clip(G*(1-mask["alpha"])+mask["color"][1]*mask["alpha"], 0, 1), G)
        B = np.where(mask["array"] == True, np.clip(B*(1-mask["alpha"])+mask["color"][2]*mask["alpha"], 0, 1), B)

    return np.stack([R, G, B], axis = 2) 

# if __name__ == "__main__":

#     # Define some paths.
#     img = {
#         "meta": {"title": "LS8 Collection 2", "collection": 2, "landsat": 8, "level": 2},
#         "qa": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LC08_L2SP_200028_20211004_20211013_02_T1/LC08_L2SP_200028_20211004_20211013_02_T1_QA_PIXEL.TIF",
#         "red": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LC08_L2SP_200028_20211004_20211013_02_T1/LC08_L2SP_200028_20211004_20211013_02_T1_SR_B4.TIF",
#         "green": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LC08_L2SP_200028_20211004_20211013_02_T1/LC08_L2SP_200028_20211004_20211013_02_T1_SR_B3.TIF",
#         "blue": r"/Users/hmcoerver/pywapor_notebooks/my_landsat_folder/LC08_L2SP_200028_20211004_20211013_02_T1/LC08_L2SP_200028_20211004_20211013_02_T1_SR_B2.TIF",
#         }

#     # Only plot [3000:4000, 4000:5000] of original image.
#     slc = [slice(3000, 4000), slice(4000, 5000)]
    
#     # Open the images as numpy arrays.
#     img_arrays = {name: open_array(path, slc = slc) for name, path in 
#                     img.items() if name != "meta"}

#     # Get collection number.
#     collection = img["meta"]["collection"]
#     ls_number = img["meta"]["landsat"]
#     level = img["meta"]["level"]

#     # Dictionary to convert between bits and flags.
#     flag_bits = get_pixel_qa_bits(collection, ls_number, level)
#     # flag_bits = get_radsat_qa_bits(collection, ls_number, level)

#     # Which flags to set to True in mask.
#     flags = ["dilated_cloud", "cloud"]
#     # flags = ["saturated_band4", "saturated_band3", "saturated_band2"]

#     # Calculate the mask.
#     my_mask = get_mask(img_arrays["qa"], flags, flag_bits)

#     # Create an RGB array.
#     mask_plot_info = {
#                         "array": my_mask,
#                         "color": (1,0,0),
#                         "alpha": 0.4,
#                         }

#     rgb = create_masked_rgb(img_arrays["red"], 
#                             img_arrays["green"], 
#                             img_arrays["blue"], mask_plot_info)

#     # Plot the image with the masks.
#     plt.figure(collection, figsize = (8,8))
#     plt.title(img["meta"]["title"])
#     plt.imshow(rgb)