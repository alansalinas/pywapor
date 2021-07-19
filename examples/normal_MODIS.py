"""
The Normal MODIS run can be used if daily data 
is required based on MODIS data. No temporal 
or spatial interpolation is used, and clouds 
in the input data are removed and not replaced. 
The run exists in two parts. Part one is to download 
and process all the required data. This will create 
the input folder. The second part will run the ETLook 
model by using this input folder. The data is stored in 
the Output folder.
"""
import os
import pandas as pd
import pywapor

# inputs
Startdate = "2019-07-06"
Enddate = "2019-07-07"
latlim = [28.5, 31.9]
lonlim = [29.2, 32.5]
output_folder = r"/Users/hmcoerver/Downloads/normal_MODIS"

dates = pd.date_range(Startdate, Enddate, freq = "D")

# Collect data for the whole daily MODIS period
pywapor.pre_et_look.main(
                        output_folder, 
                        Startdate, 
                        Enddate, 
                        latlim, 
                        lonlim,
                        LandCover = "GlobCover", 
                        Short_Downwards_Radiation = "MERRA",
                        composite = False, 
                        RAW_folder= os.path.join(output_folder, "RAW"),
                       )

# Run pyWAPOR model over daily MODIS data
for date in dates:
    pywapor.et_look_code.main(os.path.join(output_folder, "ETLook_input"), 
                                    os.path.join(output_folder, "ETLook_output"),
                                    date)