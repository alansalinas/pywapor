"""
The LS composite FRAME run can be used if an 8-daily 
composite product is required based on LS data. 
Temporal and spatial interpolation is applied, and 
clouds in the input data are removed and replaced. 
The interpolation is done as described by FRAME. The same 
steps are required as the MODIS composite FRAME run, but 
now prior to this run the LS data has to be downloaded 
and processed by the pyWAPOR.Atmospheric_Correction_LS.py 
function as described in the Normal LS run.
"""

import glob
import os
import datetime
import pyWAPOR
import watertools

# inputs
Startdate = "2019-07-03"
Enddate = "2019-07-30"
latlim = [28.7698, 29.7576]
lonlim = [30.2393, 31.2350]
output_folder = r"/content/sample_data"
Satellite_folder = r"/content/sample_data/LS_Input_Data"

# Find LS dates
os.chdir(os.path.join(Satellite_folder, "NDVI"))
files_Output = glob.glob("*.tif")
Dates = [datetime.datetime.strptime(k, '*_%Y%m%d.tif') for k in files_Output]

for Date in Dates:

    Startdate = Date.strftime("%Y-%m-%d")
    Enddate = Date.strftime("%Y-%m-%d")
    
    # Collect data for 1 day LS run
    pyWAPOR.Pre_ETLook.main(r"F:\Project_FAOJORDAN\Input_Data\Test", 
                            Startdate, 
                            Enddate, 
                            latlim,
                            lonlim, 
                            LandCover = "GlobCover", 
                            Short_Downwards_Radiation = "MERRA", 
                            RAW_folder=r"/content/sample_data/RAW",
                            composite = False, 
                            Satellite_folder = Satellite_folder)
    
    # Run pyWAPOR model for LS day
    pyWAPOR.ETLook.ETLook_code.main(os.path.join(output_folder, "ETLook_input"), 
                                    os.path.join(output_folder, "ETLook_output"), 
                                    Date)

Startdate_comp = ""
Enddate_comp = ""

# Collect data for composite LS run        
pyWAPOR.Calculate_composite_LS.main(output_folder, 
                                    Startdate_comp,
                                    Enddate_comp,
                                    latlim,
                                    lonlim,
                                    os.path.join(output_folder, "ETLook_input"), 
                                    os.path.join(output_folder, "ETLook_output"), 
                                    LandCover = "GlobCover",
                                    Short_Downwards_Radiation = "MERRA", 
                                    RAW_folder=r"/content/sample_data/RAW") 

dates_8d = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate_comp, Enddate_comp)

for date in dates_8d:
    pyWAPOR.ETLook.ETLook_code.main(os.path.join(output_folder, "ETLook_input_composite_LS"), 
                                    os.path.join(output_folder, "ETLook_output_composite_LS_FRAME"),
                                    date)