"""
The Normal LS run can be used if 8-daily data is required based 
on LS data. No temporal or spatial interpolation is used, and 
clouds in the input data are removed and not replaced. And can 
only be runned on LS overpass days.
The downloading of LS goes manually. Therefor, go first 
to https://earthexplorer.usgs.gov/ to download the required 
LS files. 
After downloading the required LS files three steps are 
required to get the WAPOR LS results.
First, use the pyWAPOR.Atmospheric_Correction_LS.py 
function to process the RAW LS files. This will produce 
a folder including the atmopheric corrected Albedo, LST, 
and NDVI of LS.
This folder needs to be defined in the Satellite_folder 
parameter in step 2. Step 3 will calculate the WAPOR results.
"""

import glob
import os
import datetime
import pyWAPOR

# inputs
latlim = [28.7698, 29.7576]
lonlim = [30.2393, 31.2350]
output_folder = r"/Users/hmcoerver/Downloads/normal_LS"
Satellite_folder = r"/Users/hmcoerver/Downloads/normal_LS/177040"

# Find LS dates
os.chdir(os.path.join(Satellite_folder, "NDVI"))
files_Output = glob.glob("*.tif")
Dates = [datetime.datetime.strptime(k.split("_")[-1], '%Y%m%d.tif') for k in files_Output]

for Date in Dates:
    
    Startdate = Date.strftime("%Y-%m-%d")
    Enddate = Date.strftime("%Y-%m-%d")
    
    # Collect date for the LS day
    pyWAPOR.Pre_ETLook.main(output_folder, 
                            Startdate, 
                            Enddate, 
                            latlim, 
                            lonlim, 
                            LandCover = "GlobCover", 
                            Short_Downwards_Radiation = "MERRA", 
                            Satellite_folder = Satellite_folder, 
                            composite = False, 
                            RAW_folder= os.path.join(output_folder, "RAW"))
    
    # Run pyWAPOR model for LS day
    pyWAPOR.ETLook.ETLook_code.main(os.path.join(output_folder, "ETLook_input"), 
                                    os.path.join(output_folder, "ETLook_output"), 
                                    Date)