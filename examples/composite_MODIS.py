"""
The MODIS composite FRAME run can be used if an 8-daily 
composite product is required based on MODIS data. 
Temporal and spatial interpolation is applied, and 
clouds in the input data are removed and replaced. 
The interpolation is done as described by FRAME. 
The relative soil moisture is used for the interpolation. 
In total 4 steps are required. The first two steps are the 
same as the Normal MODIS run. After the Normal MODIS run the 
calculated relative soil moisture is used (output of WAPOR) to 
apply the interpolation in step 3. After this interpolation the 
ETLook model is runned another time but now with 
the interpolated data.
"""
import os 
import pandas as pd
import datetime
import pywapor
import watertools

# inputs
Startdate = "2019-07-05"
Enddate = "2019-07-11"
latlim = [28.7698, 29.7576]
lonlim = [30.2393, 31.2350]
output_folder = r"/Users/hmcoerver/Downloads/composite_MODIS2"

dates_8d = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)
dates = pd.date_range(dates_8d[0], (dates_8d[-1] + datetime.timedelta(days=7)).strftime("%Y-%m-%d"), freq = "D")
Enddate = dates[-1].strftime("%Y-%m-%d")
Startdate = dates[0].strftime("%Y-%m-%d")

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
print(dates)
for date in dates:
    pywapor.et_look_code.main(os.path.join(output_folder, "ETLook_input"), 
                                    os.path.join(output_folder, "ETLook_output"),
                                    date)

pywapor.calculate_composite_modis.main(
                                       output_folder, 
                                       Startdate,
                                       Enddate,
                                       latlim,
                                       lonlim,
                                       os.path.join(output_folder, "ETLook_input"), 
                                       os.path.join(output_folder, "ETLook_output"), 
                                       LandCover = "GlobCover",
                                       Short_Downwards_Radiation = "MERRA", 
                                       RAW_folder= os.path.join(output_folder, "RAW"),
                                      )

dates_8d = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)

for date in dates_8d:
    pywapor.et_look_code.main(os.path.join(output_folder, "ETLook_input_composite_MODIS"), 
                                    os.path.join(output_folder, "ETLook_output_composite_MODIS"),
                                    date)