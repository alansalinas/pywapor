"""lue_max, landmask, z_obst_max, rs_min
"""

import numpy as np

def lulc_to_x(ds, var, convertor, out_var = None):

    new_data = ds[var] * np.nan

    for class_number in np.unique(ds[var]):
        if class_number in convertor.keys():
            new_data = new_data.where(ds[var] != class_number, 
                                        convertor[class_number])

    if "sources" in ds[var].attrs.keys():
        new_data.attrs = {"sources": ds[var].attrs["sources"]}

    if not isinstance(out_var, type(None)):
        ds[out_var] = new_data
    else: 
        ds[var] = new_data

    return ds

def wapor_to_land_mask():

    convertor = {
                    20: 1,
                    30: 1,
                    42: 1,
                    43: 1,
                    50: 3,
                    60: 1,
                    80: 2,
                    90: 1,
                    116: 1,
                    126 : 1,
                    33 : 1,
                    41 : 1,
                }

    return convertor

def wapor_to_lue_max():
    
    convertor = {
                    33: 2.66, 
                    41:	3.04,
                    42:	3.04,
                    50:	1.57,
                    60:	2.66,
                    80:	0,
    }

    return convertor

def wapor_to_rs_min():

    convertor = {
                20 : 175,
                30 : 175,
                42 : 125,
                43 : 125,
                50 : 400,
                60 : 100,
                80 : 100,
                90 : 150,
                116 : 180,
                126 : 250,
                33 : 175,
                41 : 150,
    }

    return convertor

def wapor_to_z_obst_max():

    convertor = {
                    20 : 1.0,
                    30 : 2.0,
                    42 : 1.5,
                    43 : 1.5,
                    50 : 10.0,
                    60 : 0.1,
                    80 : 0.1,
                    90 : 2.0,
                    116 : 5.0,
                    126 : 3.0,
                    33 : 8.0,
                    41 : 2.0,
    }

    return convertor

def globcover_to_land_mask():

    convertor = {
                    11: 1,    
                    14:	1,    
                    20:	1,    
                    30:	1,    
                    40:	1,    
                    50:	1,    
                    60:	1,    
                    70:	1,    
                    90:	1,    
                    100: 1,   	 
                    110: 1,   	 
                    120: 1,   	 
                    130: 1,  	 
                    140: 1,  	 
                    150: 1,  	 
                    160: 1,  	 
                    170: 1,  	 
                    180: 1,  	 
                    190: 3,  	 
                    200: 1,  	 
                    210: 2,  	 
                    220: 1,  	 
                    230: 0  	 
        }

    return convertor

def globcover_to_lue_max():
    
    convertor = {
                11: 3.04,    
                14:	3.04,    
                20:	3.04,    
                30:	3.04,    
                40:	2.45,    
                50:	2.66,    
                60:	2.66,    
                70:	2.5,    
                90:	2.5,    
                100: 2.54,   	 
                110: 2.54,   	 
                120: 2.54,   	 
                130: 2.54,  	 
                140: 2.3,  	 
                150: 1.57,  	 
                160: 2.66,  	 
                170: 1.57,  	 
                180: 1.57,  	 
                190: 1.57,  	 
                200: 1.57,  	 
                210: 0,  	 
                220: 0,  	 
                230: 0  	 
        }

    return convertor

def globcover_to_rs_min():

    convertor = {
                11: 200,    
                14:	200,    
                20:	150,    
                30:	150,    
                40:	100,    
                50:	120,    
                60:	100,    
                70:	150,    
                90:	180,    
                100: 175,   	 
                110: 150,   	 
                120: 350,   	 
                130: 175,  	 
                140: 250,  	 
                150: 150,  	 
                160: 250,  	 
                170: 200,  	 
                180: 300,  	 
                190: 100,  	 
                200: 100,  	 
                210: 100,  	 
                220: 100,  	 
                230: 0  	 
        }

    return convertor

def globcover_to_z_obst_max():

    convertor = {
            11: 4.0,    
            14:	4.0,    
            20:	2.0,    
            30:	3.5,    
            40:	0.1,    
            50:	0.6,    
            60:	1.2,    
            70:	2.0,    
            90:	5.0,    
            100: 8.0,   	 
            110: 2.0,   	 
            120: 8.0,   	 
            130: 4.0,  	 
            140: 2.0,  	 
            150: 1.0,  	 
            160: 0.3,  	 
            170: 6.0,  	 
            180: 10,  	 
            190: 10,
            200: 0.1,
            210: 0.1,
            220: 0.1, 
            230: 0,
    }

    return convertor
