# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Thu Feb 21 18:57:09 2019
"""
import os
import sys
import shutil
import datetime
import gdal
import requests
import glob
import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly

import watertools
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC

import pyWAPOR

def main(output_folder, Startdate, Enddate, latlim, lonlim, LandCover = "GlobCover", Short_Downwards_Radiation = "MSGCCP", thermal_ndvi_sharpening = False, Satellite_folder = None, composite = False, RAW_folder = None):

    ############################ Get General inputs ###############################
    
    # No composite if LS is used (Satellite folder)
    if Satellite_folder != None:
        composite = False
 
    # Define the input folders
    if composite == True:
        folders_input_RAW = os.path.join(output_folder, "RAW_composite")
        folder_input_ETLook = os.path.join(output_folder, "ETLook_input_composite")
    else:
        folders_input_RAW = os.path.join(output_folder, "RAW")
        folder_input_ETLook = os.path.join(output_folder, "ETLook_input")
 
    if RAW_folder != None:
        folders_input_RAW = RAW_folder
 
    # Create folders if not exists
    if not os.path.exists(folders_input_RAW):
        os.makedirs(folders_input_RAW)
    if not os.path.exists(folder_input_ETLook):
        os.makedirs(folder_input_ETLook)
       
    # Define the dates   
    if composite == True:
        Dates = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)
    else:
        Dates = pd.date_range(Startdate, Enddate, freq = "D")

    # Extend the days for NDVI data with +8 for both sides
    Startdate_NDVI = datetime.datetime.strptime(Startdate, "%Y-%m-%d") - datetime.timedelta(days = 8) 
    Enddate_NDVI = datetime.datetime.strptime(Enddate, "%Y-%m-%d") + datetime.timedelta(days = 8) 
    
    Startdate_NDVI_str = datetime.datetime.strftime(Startdate_NDVI, "%Y-%m-%d")
    Enddate_NDVI_str = datetime.datetime.strftime(Enddate_NDVI, "%Y-%m-%d")
    
    ######################### Download LST MODIS data #############################
    if Satellite_folder == None:
        # Download LST data
        if composite == True:
            watertools.Collect.MOD11.LST_8daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim)
            watertools.Collect.MYD11.LST_8daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim)        
            Combine_LST_composite(folders_input_RAW, Startdate, Enddate)
        else:
            watertools.Collect.MOD11.LST_daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim, angle_info = 1, time_info = 1)
            watertools.Collect.MYD11.LST_daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim, angle_info = 1, time_info = 1)   
            Combine_LST(folders_input_RAW, Startdate, Enddate)
    
    ################## Download ALBEDO and NDVI MODIS data ########################
            
        # Download NDVI and ALBEDO data
        watertools.Collect.MOD13.NDVI_16daily(folders_input_RAW, Startdate_NDVI_str, Enddate_NDVI_str, latlim, lonlim)
        watertools.Collect.MYD13.NDVI_16daily(folders_input_RAW, Startdate_NDVI_str, Enddate_NDVI_str, latlim, lonlim)
        
        if composite == True:
            
            #for Date_albedo in Dates:
            #    watertools.Collect.MCD43.Albedo_daily(folders_input_RAW, Date_albedo, Date_albedo, latlim, lonlim)
                
            watertools.Collect.MCD19.Albedo_8daily(folders_input_RAW, Startdate_NDVI_str, Enddate_NDVI_str, latlim, lonlim)
                
                
        else:     
            print("Download Albedo")
            watertools.Collect.MCD43.Albedo_daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim)
    
    ########################### Download CHIRPS data ################################    

    if composite == True:

        ######################## Download Rainfall Data ###############################
        # Download CHIRPS data
        watertools.Collect.CHIRPS.daily(folders_input_RAW, Startdate_NDVI_str, Enddate_NDVI_str, latlim, lonlim)
                
    else:     

        ######################## Download Rainfall Data ###############################
        print("Download CHIRPS")
        # Download CHIRPS data
        watertools.Collect.CHIRPS.daily(folders_input_RAW, Startdate, Enddate, latlim, lonlim)    

    ########################### Download DEM data #################################
    
    # Download DEM data
    watertools.Collect.DEM.SRTM(folders_input_RAW, latlim, lonlim)
    
    ############################ Download Landuse #################################
    if LandCover == "GlobCover":
        # Download Globcover data
        watertools.Collect.Globcover.Landuse(folders_input_RAW, latlim, lonlim)
        
    if LandCover == "WAPOR":
       
        # Download Globcover data
        watertools.Collect.WAPOR.Get_Layer(folders_input_RAW, "%s-01-01"%(Startdate.split("-")[0]), "%s-12-31"%(Enddate.split("-")[0]), latlim, lonlim, "L1_LCC_A")        

    ############################ Download Static Amplitude map #################################
    
    print("download temperature amplitudes from google drive")       
    output_folder_Tamp = os.path.join(folders_input_RAW, "GLDAS")
    if not os.path.exists(output_folder_Tamp):
        os.makedirs(output_folder_Tamp)
    
    T_amplitude_global_temp_filename = os.path.join(output_folder_Tamp, "Temp_Amplitudes_global.tif")
    if not os.path.exists(T_amplitude_global_temp_filename):
        download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", T_amplitude_global_temp_filename)

    ############### Loop over days for the dynamic data ###############################
    
    # Create the inputs of MODIS or LS for all the Dates
    for Date in Dates:
        
        try:
            # Define output folder
            folder_input_ETLook_Date = os.path.join(folder_input_ETLook, "%d%02d%02d" %(Date.year, Date.month, Date.day))
            if not os.path.exists(folder_input_ETLook_Date):
                os.makedirs(folder_input_ETLook_Date)

            folder_input_ETLook_Static = os.path.join(folder_input_ETLook, "Static")
            if not os.path.exists(folder_input_ETLook_Static):
                os.makedirs(folder_input_ETLook_Static)

            if Satellite_folder == None:
                
                # Find nearest date for NDVI 
                Startdate_year = "%d-01-01" %Date.year
                Enddate_year = "%d-12-31" %Date.year
                            
                # Create MODIS NDVI dataset
                Dates_eight_daily_year = pd.date_range(Startdate_year, Enddate_year, freq = "8D")
                
                # find nearest NDVI date
                Date_nearest = min(Dates_eight_daily_year, key=lambda Dates_eight_daily_year: abs(Dates_eight_daily_year - Date))   

                # Create NDVI files for ETLook
                
                # try MOD13 and MYD13
                NDVI_file = os.path.join(folder_input_ETLook_Date, "NDVI_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                if not os.path.exists(NDVI_file):
                    folder_RAW_file_NDVI = os.path.join(folders_input_RAW, "NDVI", "{v}13")      
                    filename_NDVI = "NDVI_{v}13Q1_-_16-daily_%d.%02d.%02d.tif" %(Date_nearest.year, Date_nearest.month, Date_nearest.day)
                    
                    if os.path.exists(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MOD")):
                        dest_ndvi = gdal.Open(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MOD"))
                        data, Geo_out, Proj_out = RC.clip_data(dest_ndvi, latlim, lonlim)
                        
                        DC.Save_as_tiff(NDVI_file, data, Geo_out, Proj_out)
                         
                    elif os.path.exists(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MYD")): 
                        
                        dest_ndvi = gdal.Open(os.path.join(folder_RAW_file_NDVI, filename_NDVI).format(v="MYD"))
                        data, Geo_out, Proj_out = RC.clip_data(dest_ndvi, latlim, lonlim)
                        
                        DC.Save_as_tiff(NDVI_file, data, Geo_out, Proj_out)
                         
                    else:
                        print("NDVI is not available for date: %d%02d%02d" %(Date.year, Date.month, Date.day))
           
            else:
                NDVI_file = os.path.join(folder_input_ETLook_Date, "NDVI_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                DEM_file = os.path.join(folders_input_RAW, "SRTM", "DEM", "DEM_SRTM_m_3s.tif")
                dest_dem = gdal.Open(DEM_file)
                geo_dem = dest_dem.GetGeoTransform()
                geo_dem = np.array(geo_dem)
                geo_dem[1] = geo_dem[1]/3
                geo_dem[5] = geo_dem[5]/3
                geo_dem = tuple(geo_dem)
                data_ex = np.ones([dest_dem.RasterYSize *3, dest_dem.RasterXSize *3]) * np.nan
                dest_dem_ex = DC.Save_as_MEM(data_ex, geo_dem, 4326)
                
                folder_NDVI_LS = os.path.join(Satellite_folder, "NDVI")
                os.chdir(folder_NDVI_LS)
                filename_NDVI = os.path.join(folder_NDVI_LS, glob.glob("*_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))[0])
                dest_rep = RC.reproject_dataset_example(filename_NDVI, dest_dem_ex, 2)
                ndvi_rep, Geo_out, Proj_Out = RC.clip_data(dest_rep, latlim, lonlim)   
                ndvi_rep[ndvi_rep==0] = -9999 
                DC.Save_as_tiff(NDVI_file, ndvi_rep, Geo_out, 4326)
        
            # Get example files
            dest_ex = gdal.Open(NDVI_file)
            geo_ex = dest_ex.GetGeoTransform()
            proj_ex = dest_ex.GetProjection()
            size_x_ex = dest_ex.RasterXSize
            size_y_ex = dest_ex.RasterYSize    

            ####################### Create lat and lon rasters ############################
            
            Lon_file = os.path.join(folder_input_ETLook_Static, "Lon.tif")     
            Lat_file = os.path.join(folder_input_ETLook_Static, "Lat.tif")            
            if not (os.path.exists(Lon_file) or os.path.exists(Lat_file)): 
                lon_deg = np.array([geo_ex[0] + np.arange(0,size_x_ex) * geo_ex[1]]*size_y_ex)
                lat_deg = np.array([geo_ex[3] + np.arange(0,size_y_ex) * geo_ex[5]]*size_x_ex).transpose()
                
                # save as tiff 
                DC.Save_as_tiff(Lon_file, lon_deg, geo_ex, proj_ex)
                DC.Save_as_tiff(Lat_file, lat_deg, geo_ex, proj_ex)
            
            else:
               dest_lon = gdal.Open(Lon_file)
               lon_deg = dest_lon.GetRasterBand(1).ReadAsArray()
               
            dlat, dlon = watertools.Functions.Area_Conversions.Area_converter.Calc_dlat_dlon(geo_ex, size_x_ex, size_y_ex)
            
            offset_GTM = int(round(lon_deg[int(lon_deg.shape[0]/2),int(lon_deg.shape[1]/2)] * 24 / 360))
        
            # Create ALBEDO files for ETLook

            if Satellite_folder == None:        
                # try MCD43
                ALBEDO_file = os.path.join(folder_input_ETLook_Date, "ALBEDO_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                if not os.path.exists(ALBEDO_file):
                    
                    if composite == True:
                        folder_RAW_file_ALBEDO = os.path.join(folders_input_RAW, "Albedo", "MCD19", "8_Daily")      
                        filename_ALBEDO = "Albedo_MCD19A3_-_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)                    
                        #folder_RAW_file_ALBEDO = os.path.join(folders_input_RAW, "Albedo", "MCD43")      
                        #filename_ALBEDO = "Albedo_MCD43A3_-_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                        
                    else:
                        folder_RAW_file_ALBEDO = os.path.join(folders_input_RAW, "Albedo", "MCD43")      
                        filename_ALBEDO = "Albedo_MCD43A3_-_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                    
                    if os.path.exists(os.path.join(folder_RAW_file_ALBEDO, filename_ALBEDO)):
                        destalbedo = RC.reproject_dataset_example(os.path.join(folder_RAW_file_ALBEDO, filename_ALBEDO), NDVI_file, method=1)
                        albedo = destalbedo.GetRasterBand(1).ReadAsArray()
                        albedo[albedo<=-0.4] = -9999                    
                        DC.Save_as_tiff(ALBEDO_file, albedo, geo_ex, proj_ex)            
                        
                    else:
                        print("ALBEDO is not available for date: %d%02d%02d" %(Date.year, Date.month, Date.day))
                        
                # Create LST files for ETLook            
                LST_file = os.path.join(folder_input_ETLook_Date, "LST_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))            
                Time_file = os.path.join(folder_input_ETLook_Date, "Time_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))       
                if not os.path.exists(LST_file):
                    if composite == True:
                        folder_RAW_file_LST = os.path.join(folders_input_RAW, "MODIS", "LST", "8_Daily")      
                        filename_LST = "LST_MCD11A2_K_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)    
                        filename_Time = "Time_MCD11A2_hour_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)                    
                    else:
                        folder_RAW_file_LST = os.path.join(folders_input_RAW, "MODIS", "LST", "Daily")      
                        filename_LST = "LST_MCD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)    
                        filename_Time = "Time_MCD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                        
                    if os.path.exists(os.path.join(folder_RAW_file_LST, filename_LST)):        
                        
                        if thermal_ndvi_sharpening == False:
                            destLST = RC.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_LST), NDVI_file, method=2)
                            LST = destLST.GetRasterBand(1).ReadAsArray()
                            LST[LST==0.0] = -9999
                            
                        else:
                            
                            # open ndvi file
                            dest_down = gdal.Open(NDVI_file)
                            NDVI = dest_down.GetRasterBand(1).ReadAsArray()
                            
                            # Create mask for thermal sharpening
                            Total_mask_thermal = np.where(NDVI>0.01, 0, 1)
                            Total_mask_thermal[Total_mask_thermal > 0] = 1
                    
                            # Open LST
                            dest_up = gdal.Open(os.path.join(folder_RAW_file_LST, filename_LST))
                            
                            # Open Thermal data
                            surface_temp_up = dest_up.GetRasterBand(1).ReadAsArray()
                            dest_lst_down = RC.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_LST), NDVI_file, method = 1)                                       
                            Surface_temp = dest_lst_down.GetRasterBand(1).ReadAsArray()
                            
                            Total_mask_thermal = np.where(Surface_temp<200, 1, Total_mask_thermal)
                            
                        	# Upscale DEM
                            Box = 7
         
                            # Upscale NDVI data
                            dest_ndvi_up = RC.reproject_dataset_example(dest_down, dest_up)
                            NDVI_up = dest_ndvi_up.GetRasterBand(1).ReadAsArray()
                    
                            # upscale the mask to coarser resolution
                            Total_mask_thermal_up = RC.resize_array_example(Total_mask_thermal, NDVI_up, method=2)
                            Total_mask_thermal_up[Total_mask_thermal_up>0]=1

                            # Remove wrong values
                            surface_temp_up[surface_temp_up<=0] = np.nan
                            NDVI_up[NDVI_up==0] = np.nan
                            surface_temp_up[surface_temp_up==1] = np.nan
                            NDVI_up[Total_mask_thermal_up==1] = np.nan
                            NDVI[Total_mask_thermal==1] = np.nan
                    
                            # Apply thermal sharpening
                            LST = Thermal_Sharpening_Linear_Forced(surface_temp_up, NDVI_up, NDVI, Box, dest_up, NDVI_file, dest_down)
                    
                            # Replace water values to original thermal  values
                            LST[np.logical_and(NDVI>-1, NDVI < 0)] = Surface_temp[np.logical_and(NDVI>-1, NDVI < 0)]
                            LST[np.isnan(LST)] = Surface_temp[np.isnan(LST)]
                            LST[LST<200] = np.nan

                        DC.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)
            
                        destTime = RC.reproject_dataset_example(os.path.join(folder_RAW_file_LST, filename_Time), NDVI_file, method=1)
                        Time = destTime.GetRasterBand(1).ReadAsArray()
                        Time[Time==0.0] = -9999
                        DC.Save_as_tiff(Time_file, Time, geo_ex, proj_ex)
                        
                    else:
                        print("LST is not available for date: %d%02d%02d" %(Date.year, Date.month, Date.day))        
                else:
                    destTime = gdal.Open(Time_file)
                    Time = destTime.GetRasterBand(1).ReadAsArray()
                    Time[Time==0.0] = -9999
            
            else:
                ALBEDO_file = os.path.join(folder_input_ETLook_Date, "ALBEDO_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
                folder_ALBEDO_LS = os.path.join(Satellite_folder, "Albedo")
                os.chdir(folder_ALBEDO_LS)
                filename_ALBEDO = os.path.join(folder_ALBEDO_LS, glob.glob("*_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))[0])
                destalbedo = RC.reproject_dataset_example(os.path.join(folder_ALBEDO_LS, filename_ALBEDO), NDVI_file, method=1)
                albedo = destalbedo.GetRasterBand(1).ReadAsArray()
                albedo[albedo<=-0.4] = -9999   
                albedo[albedo==0] = -9999                   
                DC.Save_as_tiff(ALBEDO_file, albedo, geo_ex, proj_ex)                     
          
                LST_file = os.path.join(folder_input_ETLook_Date, "LST_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))            
                Time_file = os.path.join(folder_input_ETLook_Date, "Time_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  

                folder_LST_LS = os.path.join(Satellite_folder, "LST")
                os.chdir(folder_LST_LS)
                lst_filename = glob.glob("*_LST_*_%d%02d%02d_*.tif" %(Date.year, Date.month, Date.day))[0]
                filename_LST = os.path.join(folder_LST_LS, lst_filename)
                                   

                if thermal_ndvi_sharpening == False:
                    
                    destLST = RC.reproject_dataset_example(os.path.join(folder_LST_LS, filename_LST), NDVI_file, method=1)
                    LST = destLST.GetRasterBand(1).ReadAsArray()
                    LST[LST==0.0] = -9999
                    
                else:
                    
                    # open ndvi file
                    dest_down = gdal.Open(NDVI_file)
                    NDVI = dest_down.GetRasterBand(1).ReadAsArray()
                    
                    # Create mask for thermal sharpening
                    Total_mask_thermal = np.where(NDVI>0.01, 0, 1)
                    Total_mask_thermal[Total_mask_thermal > 0] = 1
            
                    # Open LST
                    dest_start = gdal.Open(os.path.join(folder_LST_LS, filename_LST))
                    proj_up = dest_start.GetProjection()
                    geo_start = dest_start.GetGeoTransform()
                    geo_up = tuple([geo_start[0], geo_start[1] * 3, 0, geo_start[3], 0, geo_start[5] * 3])
                    Array_up = np.ones([int(np.ceil(dest_start.RasterYSize/3)), int(np.ceil(dest_start.RasterXSize/3))]) * np.nan
                    dest_up = DC.Save_as_MEM(Array_up, geo_up, proj_up)
                    
                    # Open Thermal data
                    dest_up_lst = RC.reproject_dataset_example(os.path.join(folder_LST_LS, filename_LST), dest_up, method = 4)    
                    surface_temp_up = dest_up_lst.GetRasterBand(1).ReadAsArray()
                    dest_lst_down = RC.reproject_dataset_example(os.path.join(folder_LST_LS, filename_LST), NDVI_file, method = 1)                                       
                    Surface_temp = dest_lst_down.GetRasterBand(1).ReadAsArray()
                    
                    Total_mask_thermal = np.where(Surface_temp<200, 1, Total_mask_thermal)
                    
                	# Upscale DEM
                    Box = 4
 
                    # Upscale NDVI data
                    dest_ndvi_up = RC.reproject_dataset_example(dest_down, dest_up)
                    NDVI_up = dest_ndvi_up.GetRasterBand(1).ReadAsArray()
            
                    # upscale the mask to coarser resolution
                    Total_mask_thermal_up = RC.resize_array_example(Total_mask_thermal, NDVI_up, method=2)
                    Total_mask_thermal_up[Total_mask_thermal_up>0]=1

                    # Remove wrong values
                    surface_temp_up[surface_temp_up<=0] = np.nan
                    NDVI_up[NDVI_up==0] = np.nan
                    surface_temp_up[surface_temp_up==1] = np.nan
                    NDVI_up[Total_mask_thermal_up==1] = np.nan
                    NDVI[Total_mask_thermal==1] = np.nan
            
                    # Apply thermal sharpening
                    LST = Thermal_Sharpening_Linear_Forced(surface_temp_up, NDVI_up, NDVI, Box, dest_up, NDVI_file, dest_down)
            
                    # Replace water values to original thermal  values
                    LST[np.logical_and(NDVI>-1, NDVI < 0)] = Surface_temp[np.logical_and(NDVI>-1, NDVI < 0)]
                    LST[np.isnan(LST)] = Surface_temp[np.isnan(LST)]
                    LST[LST<200] = np.nan

                DC.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)              
                
                Time = np.ones(albedo.shape)*(int(lst_filename.split("_")[-1][0:2])+(offset_GTM) + int(lst_filename.split("_")[-1][2:4])/60)
                DC.Save_as_tiff(Time_file, Time, geo_ex, proj_ex)  

            ########################## Create Time rasters ################################
            
            # calculate overall time
            dest_time = gdal.Open(Time_file)
            Time_array = dest_time.GetRasterBand(1).ReadAsArray()
            Time_array[Time_array==-9999] = np.nan
            dtime = np.nanmean(Time_array)
            if np.isnan(dtime):
                dtime = 12
            NowTime_GMT = datetime.datetime(Date.year, Date.month, Date.day, int(np.floor(dtime)), int((dtime - np.floor(dtime))*60)) - pd.DateOffset(hours = offset_GTM) 
            
            # Get DOY
            if composite == True:
                date_middle = Date + pd.DateOffset(days=4)             
                doy = int(date_middle.strftime("%j"))
 
            else:
                doy = int(Date.strftime("%j"))

            ####################### Create DEM rasters ############################
            
            # Create DEM files for ETLook            
            DEM_file = os.path.join(folder_input_ETLook_Static, "DEM.tif")                
            if not os.path.exists(DEM_file):
                folder_RAW_file_DEM = os.path.join(folders_input_RAW, "SRTM", "DEM")      
                filename_DEM = "DEM_SRTM_m_3s.tif"    
                if os.path.exists(os.path.join(folder_RAW_file_DEM, filename_DEM)):        
                    destDEM = RC.reproject_dataset_example(os.path.join(folder_RAW_file_DEM, filename_DEM), NDVI_file, method=4)
                    DEM = destDEM.GetRasterBand(1).ReadAsArray()
                    DC.Save_as_tiff(DEM_file, DEM, geo_ex, proj_ex)
                    
                else:
                    print("DEM is not available")        
        
            ##################### Calculate SLope and Aspect ##############################
            Slope_file = os.path.join(folder_input_ETLook_Static, "Slope.tif")                
            Aspect_file = os.path.join(folder_input_ETLook_Static, "Aspect.tif")   
            if not (os.path.exists(Slope_file) and os.path.exists(Aspect_file)):
                
                # open DEM
                destDEM = gdal.Open(DEM_file)
                DEM = destDEM.GetRasterBand(1).ReadAsArray()
                
                # constants
                pixel_spacing = (np.nanmean(dlon) +np.nanmean(dlat))/2
                deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
                rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree
            
                # Calculate slope
                x, y = np.gradient(DEM, pixel_spacing, pixel_spacing)
                hypotenuse_array = np.hypot(x,y)
                slope = np.arctan(hypotenuse_array) * rad2deg
            
                # calculate aspect
                aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
                aspect = 180 + aspect
        
                # Save as tiff files
                DC.Save_as_tiff(Slope_file, slope, geo_ex, proj_ex)
                DC.Save_as_tiff(Aspect_file, aspect, geo_ex, proj_ex)    
                
            ######################### Create Rainfall file ################################
            
            P_file = os.path.join(folder_input_ETLook_Date, "Precipitation_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))            
            if not os.path.exists(P_file): 
                
                folder_RAW_file_P = os.path.join(folders_input_RAW, "Precipitation", "CHIRPS", "Daily")                    
                if composite == True:
                    
                    input_format = os.path.join(folder_RAW_file_P, "P_CHIRPS.v2.0_mm-day-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif")
                    P = Calc_Composite_METEO(input_format, Date, NDVI_file, "mean")
                    P[np.isnan(P)] = 0
                    DC.Save_as_tiff(P_file, P, geo_ex, proj_ex)                   
                else:
    
                    filename_P = "P_CHIRPS.v2.0_mm-day-1_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day)
                    destP = RC.reproject_dataset_example(os.path.join(folder_RAW_file_P, filename_P), NDVI_file, method=6)
                    P = destP.GetRasterBand(1).ReadAsArray()
                    DC.Save_as_tiff(P_file, P, geo_ex, proj_ex)

            ############################# Download METEO ##################################
            
            # Define the startdates for the METEO
            if composite == True:
                Date_end = Date + pd.DateOffset(days = 7)
                StartTime = datetime.datetime(Date.year, Date.month, Date.day, 0, 0)
                EndTime = datetime.datetime(Date_end.year, Date_end.month, Date_end.day, 23, 59)               
            else:
                StartTime = datetime.datetime(Date.year, Date.month, Date.day, 0, 0)
                EndTime = datetime.datetime(Date.year, Date.month, Date.day, 23, 59)
            
            if (Date >= datetime.datetime(2000,1,1) and Date < datetime.datetime(2017,12,1)):  
            #if (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,12,1)):         
                # find nearest Meteo time            
                DateTime = pd.date_range(StartTime, EndTime, freq="H") + pd.offsets.Minute(30)
                Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - NowTime_GMT)) 
                Period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1           
    
            else:
                # find nearest Meteo time        
                DateTime = pd.date_range(StartTime, EndTime, freq="3H") + pd.offsets.Minute(90)
                Time_nearest = min(DateTime, key=lambda DateTime: abs(DateTime - NowTime_GMT)) 
                Period = np.argwhere(DateTime ==Time_nearest)[0][0] + 1
                    
         
            # Download METEO data
            #if Date < datetime.datetime(2016,1,1): 
                # Test https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero/hourly/tavg3hr_2d_asm_Nx.ascii?u10m[20:1:24][227:1:255][384:1:420]
            if Date < datetime.datetime(2000,1,1): 
                watertools.Collect.MERRA.daily(folders_input_RAW, ['u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim)
                watertools.Collect.MERRA.three_hourly(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period), int(Period+1)])
                watertools.Collect.MERRA.daily(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
         
                str_METEO = "MERRA"
                inst_name = "three_hourly"
                day_name = "daily"
                hour_steps = 3
                file_time_inst = "3-hourly"
                Periods_METEO = [int(Period), int(Period+1)]
                 
            #elif (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,12,1)):     
            elif (Date >= datetime.datetime(2000,1,1) and Date < datetime.datetime(2017,12,1)):     

                watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim)
                watertools.Collect.MERRA.hourly_MERRA2(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period), int(Period+3)])
                watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
                str_METEO = "MERRA"   
                inst_name = "hourly_MERRA2"
                day_name = "daily_MERRA2"
                hour_steps = 1
                file_time_inst = "hourly"
                Periods_METEO = [int(Period), int(Period+3)]
                
            else:
                watertools.Collect.GEOS.daily(folders_input_RAW, ['u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim)
                watertools.Collect.GEOS.three_hourly(folders_input_RAW, ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim, [int(Period), int(Period+1)])
                watertools.Collect.GEOS.daily(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
                str_METEO = "GEOS"
                inst_name = "three_hourly"
                day_name = "daily"
                hour_steps = 3
                file_time_inst = "3-hourly"
                Periods_METEO = [int(Period), int(Period+1)]
                
            # Air pressure
            pair_inst_file = os.path.join(folder_input_ETLook_Date, "Pair_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))    
            pair_inst_0_file = os.path.join(folder_input_ETLook_Date, "Pair_inst_0_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))    
            pair_24_0_file = os.path.join(folder_input_ETLook_Date, "Pair_24_0_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
                     
            if not (os.path.exists(pair_inst_file) and os.path.exists(pair_inst_0_file) and os.path.exists(pair_24_0_file)):
                folder_RAW_file_pair_inst = os.path.join(folders_input_RAW, str_METEO, "Surface_Pressure", inst_name)      
                folder_RAW_file_pair_inst_0 = os.path.join(folders_input_RAW, str_METEO, "Sea_Level_Pressure", inst_name)   
                folder_RAW_file_pair_24_0 = os.path.join(folders_input_RAW, str_METEO, "Sea_Level_Pressure", day_name) 
                HourPeriods = hour_steps * (np.array(Periods_METEO) - 1)
        
                input_format_pair_inst_MOD = os.path.join(folder_RAW_file_pair_inst, "ps_%s_kpa_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[0]))
                input_format_pair_inst_MYD = os.path.join(folder_RAW_file_pair_inst, "ps_%s_kpa_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[1]))
                    
                if composite == True:
 
                    pair_inst_MOD = Calc_Composite_METEO(input_format_pair_inst_MOD, Date, NDVI_file, "mean")
                    pair_inst_MOD[np.isnan(pair_inst_MOD)] = 0
                    
                    pair_inst_MYD = Calc_Composite_METEO(input_format_pair_inst_MYD, Date, NDVI_file, "mean")
                    pair_inst_MYD[np.isnan(pair_inst_MYD)] = 0                        
                    
                    pair_inst = np.where(Time<12, pair_inst_MOD, pair_inst_MYD)
                    
                    DC.Save_as_tiff(pair_inst_file, pair_inst, geo_ex, proj_ex)  
                    
                else:
          
                    destPairInst_MOD = RC.reproject_dataset_example(input_format_pair_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    pair_inst_MOD = destPairInst_MOD.GetRasterBand(1).ReadAsArray()

                    destPairInst_MYD = RC.reproject_dataset_example(input_format_pair_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    pair_inst_MYD = destPairInst_MYD.GetRasterBand(1).ReadAsArray()

                    pair_inst = np.where(Time<12, pair_inst_MOD, pair_inst_MYD)

                    DC.Save_as_tiff(pair_inst_file, pair_inst, geo_ex, proj_ex)
 
        
                input_format_pair_inst_sea_MOD = os.path.join(folder_RAW_file_pair_inst_0, "slp_%s_kpa_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[0]))
                input_format_pair_inst_sea_MYD = os.path.join(folder_RAW_file_pair_inst_0, "slp_%s_kpa_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[1]))

                if composite == True:
                    
                    pair_inst_sea_MOD = Calc_Composite_METEO(input_format_pair_inst_sea_MOD, Date, NDVI_file, "mean")
                    pair_inst_sea_MOD[np.isnan(pair_inst_sea_MOD)] = 0
                    
                    pair_inst_sea_MYD = Calc_Composite_METEO(input_format_pair_inst_sea_MYD, Date, NDVI_file, "mean")
                    pair_inst_sea_MYD[np.isnan(pair_inst_sea_MYD)] = 0                        
                    
                    pair_inst_sea = np.where(Time<12, pair_inst_sea_MOD, pair_inst_sea_MYD)
                    
                    DC.Save_as_tiff(pair_inst_0_file, pair_inst_sea, geo_ex, proj_ex)  
                else:    
                    
                    destPairInstSea_MOD = RC.reproject_dataset_example(input_format_pair_inst_sea_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    pair_inst_sea_MOD = destPairInstSea_MOD.GetRasterBand(1).ReadAsArray()

                    destPairInstSea_MYD = RC.reproject_dataset_example(input_format_pair_inst_sea_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    pair_inst_sea_MYD = destPairInstSea_MYD.GetRasterBand(1).ReadAsArray()

                    Pair_inst_sea = np.where(Time<12, pair_inst_sea_MOD, pair_inst_sea_MYD)

                    DC.Save_as_tiff(pair_inst_0_file, Pair_inst_sea, geo_ex, proj_ex)                        
  
       
                input_format_pair_24_sea = os.path.join(folder_RAW_file_pair_24_0, "slp_%s_kpa_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))                   
                if composite == True:
                    
                    pair_24_sea = Calc_Composite_METEO(input_format_pair_24_sea, Date, NDVI_file, "mean")
                    pair_24_sea[np.isnan(pair_24_sea)] = 0                       
                    
                    DC.Save_as_tiff(pair_24_0_file, pair_24_sea, geo_ex, proj_ex)
                    
                else:                
                    
                    destPair24Sea = RC.reproject_dataset_example(input_format_pair_24_sea.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    Pair_24_sea = destPair24Sea.GetRasterBand(1).ReadAsArray()
                    DC.Save_as_tiff(pair_24_0_file, Pair_24_sea, geo_ex, proj_ex)
                    
                    
            # Specific Humidity
            qv_inst_file = os.path.join(folder_input_ETLook_Date, "qv_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))    
            qv_24_file = os.path.join(folder_input_ETLook_Date, "qv_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
                     
            if not (os.path.exists(qv_inst_file) and os.path.exists(qv_24_file)):
                folder_RAW_file_qv_inst = os.path.join(folders_input_RAW, str_METEO, "Specific_Humidity", inst_name)      
                folder_RAW_file_qv_24 = os.path.join(folders_input_RAW, str_METEO, "Specific_Humidity", day_name) 
                HourPeriods = hour_steps * (np.array(Periods_METEO) - 1)
                if str_METEO == "MERRA":
                    para = "q2m"
                else:
                    para = "qv2m"
 
                input_format_qv_inst_MOD = os.path.join(folder_RAW_file_qv_inst, "%s_%s_kg-kg-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(para, str_METEO, file_time_inst, HourPeriods[0]))
                input_format_qv_inst_MYD = os.path.join(folder_RAW_file_qv_inst, "%s_%s_kg-kg-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(para, str_METEO, file_time_inst, HourPeriods[1]))
                input_format_qv_24 = os.path.join(folder_RAW_file_qv_24, "%s_%s_kg-kg-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(para, str_METEO))

                if composite == True:
   
                    qv_inst_MOD = Calc_Composite_METEO(input_format_qv_inst_MOD, Date, NDVI_file, "mean")
                    qv_inst_MOD[np.isnan(qv_inst_MOD)] = 0
                    
                    qv_inst_MYD = Calc_Composite_METEO(input_format_qv_inst_MYD, Date, NDVI_file, "mean")
                    qv_inst_MYD[np.isnan(qv_inst_MYD)] = 0                        
                    
                    qv_inst = np.where(Time<12, qv_inst_MOD, qv_inst_MYD)
                    
                    DC.Save_as_tiff(qv_inst_file, qv_inst, geo_ex, proj_ex)  
             
                else:
   
                    destqvInst_MOD = RC.reproject_dataset_example(input_format_qv_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    qv_inst_MOD = destqvInst_MOD.GetRasterBand(1).ReadAsArray()
                    
                    destqvInst_MYD = RC.reproject_dataset_example(input_format_qv_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    qv_inst_MYD = destqvInst_MYD.GetRasterBand(1).ReadAsArray()                    
 
                    qv_inst = np.where(Time<12, qv_inst_MOD, qv_inst_MYD)
                   
                    DC.Save_as_tiff(qv_inst_file, qv_inst, geo_ex, proj_ex)

                if composite == True:
                    
                    qv_24 = Calc_Composite_METEO(input_format_qv_24, Date, NDVI_file, "mean")
                    qv_24[np.isnan(qv_24)] = 0                       
                    
                    DC.Save_as_tiff(qv_24_file, qv_24, geo_ex, proj_ex)
                    
                else:

                    destqv24 = RC.reproject_dataset_example(input_format_qv_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    qv_24 = destqv24.GetRasterBand(1).ReadAsArray()
                    DC.Save_as_tiff(qv_24_file, qv_24, geo_ex, proj_ex)

            # Air temperature
            Tair_inst_file = os.path.join(folder_input_ETLook_Date, "tair_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))    
            Tair_24_file = os.path.join(folder_input_ETLook_Date, "tair_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
            Tair_max_24_file = os.path.join(folder_input_ETLook_Date, "tair_max_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day)) 
            Tair_min_24_file = os.path.join(folder_input_ETLook_Date, "tair_min_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day)) 
                     
            if not (os.path.exists(Tair_inst_file) and os.path.exists(Tair_24_file) and os.path.exists(Tair_max_24_file) and os.path.exists(Tair_min_24_file)):
                folder_RAW_file_tair_inst = os.path.join(folders_input_RAW, str_METEO, "Air_Temperature", inst_name)      
                folder_RAW_file_tair_24 = os.path.join(folders_input_RAW, str_METEO, "Air_Temperature", day_name) 
                HourPeriods = hour_steps * (np.array(Periods_METEO) - 1)
                
                input_format_tair_inst_MOD = os.path.join(folder_RAW_file_tair_inst, "t2m_%s_K_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[0]))
                input_format_tair_inst_MYD = os.path.join(folder_RAW_file_tair_inst, "t2m_%s_K_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[1]))
                input_format_tair_24 = os.path.join(folder_RAW_file_tair_24, "t2m_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
                input_format_tair_min_24 = os.path.join(folder_RAW_file_tair_24,  "min", "t2mmin_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
                input_format_tair_max_24 = os.path.join(folder_RAW_file_tair_24,  "max", "t2mmax_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))

                if composite == True:

                    tair_inst_MOD = Calc_Composite_METEO(input_format_tair_inst_MOD, Date, NDVI_file, "max", DEM_file = DEM_file, lapse = -0.008)
                    tair_inst_MOD[np.isnan(tair_inst_MOD)] = 0
                    
                    tair_inst_MYD = Calc_Composite_METEO(input_format_tair_inst_MYD, Date, NDVI_file, "max", DEM_file = DEM_file, lapse = -0.008)
                    tair_inst_MYD[np.isnan(tair_inst_MYD)] = 0                        
                    
                    tair_inst = np.where(Time<12, tair_inst_MOD, tair_inst_MYD)
                    if np.nanmax(tair_inst>270):
                        tair_inst = tair_inst -273.15
                    
                    DC.Save_as_tiff(Tair_inst_file, tair_inst, geo_ex, proj_ex)                      

                else:
                    
                    #destTairInst_MOD = RC.reproject_dataset_example(input_format_tair_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    tair_inst_MOD = lapse_rate_temp(input_format_tair_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), DEM_file, lapse = -0.008)
                    
                    #destTairInst_MYD = RC.reproject_dataset_example(input_format_tair_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    tair_inst_MYD = lapse_rate_temp(input_format_tair_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), DEM_file, lapse = -0.008)                 
 
                    tair_inst = np.where(Time<12, tair_inst_MOD, tair_inst_MYD)
                    if np.nanmax(tair_inst>270):
                        tair_inst = tair_inst -273.15
                        
                    DC.Save_as_tiff(Tair_inst_file, tair_inst, geo_ex, proj_ex)  
        
                if composite == True:
                    
                    tair_24 = Calc_Composite_METEO(input_format_tair_24, Date, NDVI_file, "max", DEM_file = DEM_file, lapse = -0.006)
                    tair_24[np.isnan(tair_24)] = 0                       
                    if np.nanmax(tair_24>270):
                        tair_24 = tair_24 -273.15
                        
                    DC.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)
                    
                else:

                    #desttair24 = RC.reproject_dataset_example(input_format_tair_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    tair_24 = lapse_rate_temp(input_format_tair_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), DEM_file, lapse = -0.006)
                    if np.nanmax(tair_24>270):
                        tair_24 = tair_24 -273.15

                    DC.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)
                    

                if composite == True:
                    
                    tair_min_24 = Calc_Composite_METEO(input_format_tair_min_24, Date, NDVI_file, "max", DEM_file = DEM_file, lapse = 0.0)
                    tair_min_24[np.isnan(tair_min_24)] = 0                       
                    if np.nanmax(tair_min_24>270):
                        tair_min_24 = tair_min_24 -273.15
                                           
                    DC.Save_as_tiff(Tair_min_24_file, tair_min_24, geo_ex, proj_ex)                   
                    
                else:

                    #desttairmin24 = RC.reproject_dataset_example(input_format_tair_min_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    tair_min_24 = lapse_rate_temp(input_format_tair_min_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), DEM_file, lapse = 0.0)
                    if np.nanmax(tair_min_24>270):
                        tair_min_24 = tair_min_24 -273.15                    
                    DC.Save_as_tiff(Tair_min_24_file, tair_min_24, geo_ex, proj_ex)                    
                    
                    
                if composite == True:
                    
                    tair_max_24 = Calc_Composite_METEO(input_format_tair_max_24, Date, NDVI_file, "max", DEM_file = DEM_file, lapse = 0.0)
                    tair_max_24[np.isnan(tair_max_24)] = 0                       
                    if np.nanmax(tair_max_24>270):
                        tair_max_24 = tair_max_24 -273.15                       
                    DC.Save_as_tiff(Tair_max_24_file, tair_max_24, geo_ex, proj_ex)                   
                    
                else:

                    #desttairmax24 = RC.reproject_dataset_example(input_format_tair_max_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    tair_max_24 = lapse_rate_temp(input_format_tair_max_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), DEM_file, lapse = 0.0)
                    if np.nanmax(tair_max_24>270):
                        tair_max_24 = tair_max_24 -273.15                           
                    DC.Save_as_tiff(Tair_max_24_file, tair_max_24, geo_ex, proj_ex)                    
                    
        
            # Wind Speed
            wind_inst_file = os.path.join(folder_input_ETLook_Date, "wind_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))    
            wind_24_file = os.path.join(folder_input_ETLook_Date, "wind_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
                     
            if not (os.path.exists(wind_inst_file) and os.path.exists(wind_24_file)):
                
                HourPeriods = hour_steps * (np.array(Periods_METEO) - 1)
                
                folder_RAW_file_u2_inst = os.path.join(folders_input_RAW, str_METEO, "Eastward_Wind", inst_name)      
                folder_RAW_file_u2_24 = os.path.join(folders_input_RAW, str_METEO, "Eastward_Wind", day_name) 
                folder_RAW_file_v2_inst = os.path.join(folders_input_RAW, str_METEO, "Northward_Wind", inst_name)      
                folder_RAW_file_v2_24 = os.path.join(folders_input_RAW, str_METEO, "Northward_Wind", day_name) 

                input_format_u2_inst_MOD = os.path.join(folder_RAW_file_u2_inst, "u2m_%s_m-s-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[0]))
                input_format_u2_inst_MYD = os.path.join(folder_RAW_file_u2_inst, "u2m_%s_m-s-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[1]))
                input_format_v2_inst_MOD = os.path.join(folder_RAW_file_v2_inst, "v2m_%s_m-s-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[0]))
                input_format_v2_inst_MYD = os.path.join(folder_RAW_file_v2_inst, "v2m_%s_m-s-1_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(str_METEO, file_time_inst, HourPeriods[1]))
                input_format_u2_24 = os.path.join(folder_RAW_file_u2_24, "u2m_%s_m-s-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
                input_format_v2_24 = os.path.join(folder_RAW_file_v2_24, "v2m_%s_m-s-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))


                if composite == True:
                    u2_inst_MOD = Calc_Composite_METEO(input_format_u2_inst_MOD, Date, NDVI_file, "mean")
                    u2_inst_MOD[np.isnan(u2_inst_MOD)] = 0
                    
                    u2_inst_MYD = Calc_Composite_METEO(input_format_u2_inst_MYD, Date, NDVI_file, "mean")
                    u2_inst_MYD[np.isnan(u2_inst_MYD)] = 0                        
                    
                    u2_inst = np.where(Time<12, u2_inst_MOD, u2_inst_MYD)
                    
                    v2_inst_MOD = Calc_Composite_METEO(input_format_v2_inst_MOD, Date, NDVI_file, "mean")
                    v2_inst_MOD[np.isnan(v2_inst_MOD)] = 0
                    
                    v2_inst_MYD = Calc_Composite_METEO(input_format_v2_inst_MYD, Date, NDVI_file, "mean")
                    v2_inst_MYD[np.isnan(v2_inst_MYD)] = 0                        
                    
                    v2_inst = np.where(Time<12, v2_inst_MOD, v2_inst_MYD)
                    
                    wind_inst = np.sqrt(u2_inst**2 + v2_inst **2)                   

                    DC.Save_as_tiff(wind_inst_file, wind_inst, geo_ex, proj_ex)                   
                    
                else:

                    destu2inst_MOD = RC.reproject_dataset_example(input_format_u2_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    destu2inst_MYD = RC.reproject_dataset_example(input_format_u2_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    destv2inst_MOD = RC.reproject_dataset_example(input_format_v2_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    destv2inst_MYD = RC.reproject_dataset_example(input_format_v2_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)

                    u2_inst_MOD = destu2inst_MOD.GetRasterBand(1).ReadAsArray()
                    u2_inst_MYD = destu2inst_MYD.GetRasterBand(1).ReadAsArray()
                    v2_inst_MOD = destv2inst_MOD.GetRasterBand(1).ReadAsArray()
                    v2_inst_MYD = destv2inst_MYD.GetRasterBand(1).ReadAsArray()
                    
                    u2_inst = np.where(Time<12, u2_inst_MOD, u2_inst_MYD)                    
                    v2_inst = np.where(Time<12, v2_inst_MOD, v2_inst_MYD)
                    
                    wind_inst = np.sqrt(u2_inst**2 + v2_inst **2)                       
                    
                    DC.Save_as_tiff(wind_inst_file, wind_inst, geo_ex, proj_ex)

                if composite == True:                   
                    
                    u2_24= Calc_Composite_METEO(input_format_u2_24, Date, NDVI_file, "mean")
                    u2_24[np.isnan(u2_24)] = 0                       
 
                    v2_24= Calc_Composite_METEO(input_format_v2_24, Date, NDVI_file, "mean")
                    v2_24[np.isnan(v2_24)] = 0    
                    
                    wind_24 = np.sqrt(u2_24**2 + v2_24 **2) 
                     
                    DC.Save_as_tiff(wind_24_file, wind_24, geo_ex, proj_ex)                      
                    
                else:
                     
                    destu224 = RC.reproject_dataset_example(input_format_u2_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    u2_24 = destu224.GetRasterBand(1).ReadAsArray()
                    
                    destv224 = RC.reproject_dataset_example(input_format_v2_24.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    v2_24 = destv224.GetRasterBand(1).ReadAsArray()                    
 
                    wind_24 = np.sqrt(u2_24**2 + v2_24 **2) 
                   
                    DC.Save_as_tiff(wind_24_file, wind_24, geo_ex, proj_ex)
                    
            wv_inst_file = os.path.join(folder_input_ETLook_Date, "wv_inst_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))           
            if not os.path.exists(wv_inst_file):  
                
                HourPeriods = hour_steps * (np.array(Periods_METEO) - 1)
                if str_METEO == "MERRA":
                    para = "tpw"
                else:
                    para = "tqv"
                                           
                folder_RAW_file_wv_inst = os.path.join(folders_input_RAW, str_METEO, "Total_Precipitable_Water_Vapor", inst_name) 
                input_format_wv_inst_MOD = os.path.join(folder_RAW_file_wv_inst, "%s_%s_mm_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(para, str_METEO, file_time_inst, HourPeriods[0]))
                input_format_wv_inst_MYD = os.path.join(folder_RAW_file_wv_inst, "%s_%s_mm_%s_{yyyy}.{mm:02d}.{dd:02d}_H%02d.M00.tif" %(para, str_METEO, file_time_inst, HourPeriods[1]))

                if composite == True:
                    
                    wv_inst_MOD = Calc_Composite_METEO(input_format_wv_inst_MOD, Date, NDVI_file, "mean")
                    wv_inst_MOD[np.isnan(wv_inst_MOD)] = 0
                    
                    wv_inst_MYD = Calc_Composite_METEO(input_format_wv_inst_MYD, Date, NDVI_file, "mean")
                    wv_inst_MYD[np.isnan(wv_inst_MYD)] = 0                        
                    
                    wv_inst = np.where(Time<12, wv_inst_MOD, wv_inst_MYD)
                    
                    DC.Save_as_tiff(wv_inst_file, wv_inst, geo_ex, proj_ex)                      

                else:    

                    destwvInst_MOD = RC.reproject_dataset_example(input_format_wv_inst_MOD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    wv_inst_MOD = destwvInst_MOD.GetRasterBand(1).ReadAsArray()
                    
                    destwvInst_MYD = RC.reproject_dataset_example(input_format_wv_inst_MYD.format(yyyy=Date.year, mm=Date.month, dd=Date.day), NDVI_file, method=6)
                    wv_inst_MYD = destwvInst_MYD.GetRasterBand(1).ReadAsArray()     
                    
                    wv_inst = np.where(Time<12, wv_inst_MOD, wv_inst_MYD)
                    
                    DC.Save_as_tiff(wv_inst_file, wv_inst, geo_ex, proj_ex)                    
                    

            ##################### Calculate Landmask ##############################
                      
            LM_file = os.path.join(folder_input_ETLook_Static, "LandMask_%s.tif" %Date.year)   
            Bulk_file = os.path.join(folder_input_ETLook_Static, "Bulk_Stomatal_resistance_%s.tif" %Date.year)   
            MaxObs_file = os.path.join(folder_input_ETLook_Static, "Maximum_Obstacle_Height_%s.tif" %Date.year)  
            LUEmax_file = os.path.join(folder_input_ETLook_Static, "LUEmax_%s.tif" %Date.year)  
            if not (os.path.exists(LM_file) and os.path.exists(Bulk_file) and os.path.exists(MaxObs_file) and os.path.exists(LUEmax_file)):
                
                if LandCover == "GlobCover":
                    folder_RAW_file_LC = os.path.join(folders_input_RAW, "GlobCover", "Landuse")      
                    filename_LC = "LC_GLOBCOVER_V2.3.tif"    
                if LandCover == "WAPOR":
                    folder_RAW_file_LC = os.path.join(folders_input_RAW, "L1_LCC_A")      
                    filename_LC = "L1_LCC_A_WAPOR_YEAR_%s.01.01.tif" %(Date.year)                      
                    
                if os.path.exists(os.path.join(folder_RAW_file_LC, filename_LC)):        
                    destLC = RC.reproject_dataset_example(os.path.join(folder_RAW_file_LC, filename_LC), NDVI_file, method=1)
                    LC = destLC.GetRasterBand(1).ReadAsArray()
                    LC[np.isnan(LC)] = -9999
                    
                    # import list with numbers to convert globcover into other maps
                    import pyWAPOR.Functions.LandCover_Converter as LCC
                    
                    if LandCover == "GlobCover":
                        # Get conversion between globcover and landmask
                        LU_LM_Classes = LCC.Globcover_LM()
                        LU_Bulk_Classes = LCC.Globcover_Bulk()
                        LU_MaxObs_Classes = LCC.Globcover_MaxObs()
                        LU_LUEmax_Classes = LCC.Globcover_LUEmax()
                        
                    if LandCover == "WAPOR":
                        # Get conversion between globcover and landmask
                        LU_LM_Classes = LCC.WAPOR_LM()
                        LU_Bulk_Classes = LCC.WAPOR_Bulk()
                        LU_MaxObs_Classes = LCC.WAPOR_MaxObs()    
                        LU_LUEmax_Classes = LCC.WAPOR_LUEmax()
                        
                    # Create Array for LandMask
                    LM = np.ones([size_y_ex, size_x_ex]) * np.nan            
                    Bulk = np.ones([size_y_ex, size_x_ex]) * np.nan 
                    MaxObs = np.ones([size_y_ex, size_x_ex]) * np.nan 
                    LUEmax = np.ones([size_y_ex, size_x_ex]) * np.nan 
                     
                    # Create LandMask 
                    for LU_LM_Class in LU_LM_Classes.keys():
                        Value_LM = LU_LM_Classes[LU_LM_Class]
                        Value_Bulk = LU_Bulk_Classes[LU_LM_Class]
                        Value_MaxObs = LU_MaxObs_Classes[LU_LM_Class]                
                        Value_LUEmax = LU_LUEmax_Classes[LU_LM_Class]   
                        LM[LC == LU_LM_Class] = Value_LM
                        Bulk[LC == LU_LM_Class] = Value_Bulk
                        MaxObs[LC  == LU_LM_Class] = Value_MaxObs
                        LUEmax[LC  == LU_LM_Class] = Value_LUEmax
                        
                    # Save as tiff files               
                    DC.Save_as_tiff(LM_file, LM, geo_ex, proj_ex)
                    DC.Save_as_tiff(Bulk_file, Bulk, geo_ex, proj_ex)
                    DC.Save_as_tiff(MaxObs_file, MaxObs, geo_ex, proj_ex)
                    DC.Save_as_tiff(LUEmax_file, LUEmax, geo_ex, proj_ex)
                    
                else:
                    print("LandCover is not available")       
            
            ########################### Download amplitude ################################
                
            # yearly amplitude temperature air
            Tair_amp_file = os.path.join(folder_input_ETLook_Static, "Tair_amp_%d.tif" %(Date.year))     
            output_folder_Tamp = os.path.join(folders_input_RAW, "GLDAS")       
            
            if not os.path.exists(Tair_amp_file):
                T_amplitude_global_temp_filename = os.path.join(output_folder_Tamp, "Temp_Amplitudes_global.tif")
                if os.path.exists(T_amplitude_global_temp_filename):        
                    desttairamp = RC.reproject_dataset_example(T_amplitude_global_temp_filename, NDVI_file, method=6)
                    tair_amp = desttairamp.GetRasterBand(1).ReadAsArray()
                    DC.Save_as_tiff(Tair_amp_file, tair_amp, geo_ex, proj_ex)
                    
                else:
                    print("Yearly Tair amplitude is not available")    
        
            ######################## Download Transmissivity ##############################
            
            # Download MSGCPP data
            #if Date < datetime.datetime(2016,1,1): 
            if Date < datetime.datetime(2002,1,1): 
                watertools.Collect.MERRA.daily(folders_input_RAW, ['swgnet'],StartTime, EndTime, latlim, lonlim)
                str_TRANS = "MERRA"
                day_name = "daily"

            elif (Date >= datetime.datetime(2002,1,1) and Date < datetime.datetime(2017,1,1)):        
            #elif (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,1,1)):
                watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['swgnet'],StartTime, EndTime, latlim, lonlim)
                str_TRANS = "MERRA"    
                day_name = "daily_MERRA2"
                
            else:
                if Short_Downwards_Radiation ==  "MSGCCP":
                    watertools.Collect.MSGCPP.SDS(folders_input_RAW, StartTime, EndTime, latlim, lonlim)
                    str_TRANS = "MSGCPP"
                    day_name = "daily"
                if Short_Downwards_Radiation == "MERRA":
                    watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['swgnet'],StartTime, EndTime, latlim, lonlim)
                    str_TRANS = "MERRA"    
                    day_name = "daily_MERRA2"            
                    
            # yearly amplitude temperature air
            Trans_file = os.path.join(folder_input_ETLook_Date, "Trans_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))   
                     
            if not os.path.exists(Trans_file):
                
                # Calculate the extraterrestrial daily radiation
                destLat = gdal.Open(Lat_file)
                lat = destLat.GetRasterBand(1).ReadAsArray()
                
                Gsc = 1367        # Solar constant (W / m2)
                deg2rad = np.pi / 180.0
                # Computation of Hour Angle (HRA = w)
                B = 360./365 * (doy-81)           # (degrees)
                # Computation of cos(theta), where theta is the solar incidence angle
                # relative to the normal to the land surface
                delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)
                
                phi = lat * deg2rad                                     # latitude of the pixel (radians)
                
                dr = 1 + 0.033 * np.cos(doy*2*np.pi/365)
                
                # Daily 24 hr radiation - For flat terrain only !
                ws_angle = np.arccos(-np.tan(phi)*np.tan(delta))   # Sunset hour angle ws   
                
                # Extraterrestrial daily radiation, Ra (W/m2):
                Ra24_flat = (Gsc/np.pi * dr * (ws_angle * np.sin(phi) * np.sin(delta) +
                                np.cos(phi) * np.cos(delta) * np.sin(ws_angle)))
        
                if composite == True:
                    
                    Dates_End_Rads = Date + datetime.timedelta(days = 7) 
                    Dates_rads = pd.date_range(Date, Dates_End_Rads)
                    i = 0
                    for Date_rad in Dates_rads:
                        if str_TRANS == "MERRA":
                            folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "Surface_Net_Downward_Shortwave_Flux", day_name) 
                            filename_trans = "swgnet_MERRA_W-m-2_daily_%d.%02d.%02d.tif"  %(Date.year, Date.month, Date.day)            
                        if str_TRANS == "MSGCPP":
                            folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "SDS", "15min") 
                            filename_trans = "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H{hour}.M{minutes}.tif"  %(Date.year, Date.month, Date.day)
                            
                            os.chdir(folder_RAW_file_trans)
                            files = glob.glob(filename_trans.format(hour = "*", minutes = "*"))
                            i = 0
                            
                            # Open all the 15 minutes files
                            for file in files:
                                file_in = os.path.join(folder_RAW_file_trans, file)
                                destswgone = gdal.Open(file_in)
                                try:
                                    swgnet_one = destswgone.GetRasterBand(1).ReadAsArray()
                                    swgnet_one[swgnet_one<0] = 0                 
                                    if not "geo_trans" in locals():
                                        swgnet = np.ones([destswgone.RasterYSize, destswgone.RasterXSize, len(files)]) * np.nan
                                        geo_trans = destswgone.GetGeoTransform()
                                        proj_trans = destswgone.GetProjection()
                                    swgnet[:,:,i] = swgnet_one
                                except:
                                    pass
                                i+=1 
                            
                            # Calculate the daily mean     
                            swgnet_mean = np.nanmean(swgnet, 2)
                            dest_swgnet_mean = DC.Save_as_MEM(swgnet_mean,geo_trans, proj_trans)
                            destswgnet = RC.reproject_dataset_example(dest_swgnet_mean, NDVI_file, method=6)
                            del geo_trans
                            
                        else:
                            destswgnet = RC.reproject_dataset_example(os.path.join(folder_RAW_file_trans, filename_trans), NDVI_file, method=6)
                             
                        swgnet_one = destswgnet.GetRasterBand(1).ReadAsArray()
                        if i == 0:
                            swgnet = np.ones([8, destswgnet.RasterYSize, destswgnet.RasterXSize])
                        swgnet[i, :, :] = swgnet_one
                        i+=1
                        
                    trans = np.nanmean(swgnet,axis = 0) / Ra24_flat
                    DC.Save_as_tiff(Trans_file, trans, geo_ex, proj_ex)
                                    
           
                else:
                    
                    if str_TRANS == "MERRA":
                        folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "Surface_Net_Downward_Shortwave_Flux", day_name) 
                        filename_trans = "swgnet_MERRA_W-m-2_daily_%d.%02d.%02d.tif"  %(Date.year, Date.month, Date.day)            
                    if str_TRANS == "MSGCPP":
                        folder_RAW_file_trans = os.path.join(folders_input_RAW, str_TRANS, "SDS", "15min") 
                        filename_trans = "SDS_MSGCPP_W-m-2_15min_%d.%02d.%02d_H{hour}.M{minutes}.tif"  %(Date.year, Date.month, Date.day)
                        
                        os.chdir(folder_RAW_file_trans)
                        files = glob.glob(filename_trans.format(hour = "*", minutes = "*"))
                        i = 0
                        
                        # Open all the 15 minutes files
                        for file in files:
                            file_in = os.path.join(folder_RAW_file_trans, file)
                            destswgone = gdal.Open(file_in)
                            try:
                                swgnet_one = destswgone.GetRasterBand(1).ReadAsArray()
                                swgnet_one[swgnet_one<0] = 0                 
                                if not "geo_trans" in locals():
                                    swgnet = np.ones([destswgone.RasterYSize, destswgone.RasterXSize, len(files)]) * np.nan
                                    geo_trans = destswgone.GetGeoTransform()
                                    proj_trans = destswgone.GetProjection()
                                swgnet[:,:,i] = swgnet_one
                            except:
                                pass
                            i+=1 
                        
                        # Calculate the daily mean     
                        swgnet_mean = np.nanmean(swgnet, 2)
                        dest_swgnet_mean = DC.Save_as_MEM(swgnet_mean,geo_trans, proj_trans)
                        destswgnet = RC.reproject_dataset_example(dest_swgnet_mean, NDVI_file, method=6)
                        del geo_trans
                        
                    else:
                        destswgnet = RC.reproject_dataset_example(os.path.join(folder_RAW_file_trans, filename_trans), NDVI_file, method=6)
                         
                    swgnet = destswgnet.GetRasterBand(1).ReadAsArray()
                    trans = swgnet / Ra24_flat
                    DC.Save_as_tiff(Trans_file, trans, geo_ex, proj_ex)
                
        except Exception as e:
            print("No ETLook input dataset for %s" %Date)
            print(e)  
            
    return()    
    
def lapse_rate_temp(tair_file, dem_file, lapse):
        
    destT_down = RC.reproject_dataset_example(tair_file, dem_file, 2)
    destDEM_up = RC.reproject_dataset_example(dem_file, tair_file, 4)
    destDEM_down = gdal.Open(dem_file)
    destDEM_up_down = RC.reproject_dataset_example(destDEM_up, dem_file, 2)
    
    # Open Arrays
    T = destT_down.GetRasterBand(1).ReadAsArray()
    DEM_down = destDEM_down.GetRasterBand(1).ReadAsArray()
    DEM_up_ave = destDEM_up_down.GetRasterBand(1).ReadAsArray()
    
    # correct wrong values
    DEM_down[DEM_down<=0]=0
    DEM_up_ave[DEM_up_ave<=0]=0
    
    # 
    Tdown = pyWAPOR.ETLook.meteo.disaggregate_air_temperature(T, DEM_down, DEM_up_ave, lapse)

    return(Tdown)
    
def Combine_LST(folders_input_RAW, Startdate, Enddate):
    
    Dates = pd.date_range(Startdate, Enddate, freq = "D")
    
    output_folder_end = os.path.join(folders_input_RAW, "MODIS", "LST", "Daily")
    
    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)
    
    for Date in Dates:
  
        LST_file = os.path.join(output_folder_end, "LST_MCD11A1_K_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        Time_file = os.path.join(output_folder_end, "Time_MCD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))    

        if not (os.path.exists(Time_file) or os.path.exists(LST_file)):
            filename_angle_mod = os.path.join(folders_input_RAW, "LST", "MOD11", "Daily", "Angle_MOD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_angle_myd = os.path.join(folders_input_RAW, "LST", "MYD11", "Daily", "Angle_MYD11A1_degrees_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_mod = os.path.join(folders_input_RAW, "LST", "MOD11", "Daily", "Time_MOD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            filename_time_myd = os.path.join(folders_input_RAW, "LST", "MYD11", "Daily", "Time_MYD11A1_hour_daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            format_lst_mod = os.path.join(folders_input_RAW, "LST", "MOD11", "Daily", "LST_MOD11A1_K_daily_%d.%02d.%02d.{time}.tif" %(Date.year, Date.month, Date.day))
            format_lst_myd = os.path.join(folders_input_RAW, "LST", "MYD11", "Daily", "LST_MYD11A1_K_daily_%d.%02d.%02d.{time}.tif" %(Date.year, Date.month, Date.day))
            os.chdir(os.path.join(folders_input_RAW, "LST", "MOD11", "Daily"))
            filename_lst_mod = glob.glob(format_lst_mod.format(time = "*"))[0]
            os.chdir(os.path.join(folders_input_RAW, "LST", "MYD11", "Daily"))           
            filename_lst_myd = glob.glob(format_lst_myd.format(time = "*"))[0]
                
            dest_angle_mod = gdal.Open(filename_angle_mod)   
            dest_angle_myd = gdal.Open(filename_angle_myd)        
            dest_time_mod = gdal.Open(filename_time_mod)       
            dest_time_myd = gdal.Open(filename_time_myd)  
            dest_lst_mod = gdal.Open(filename_lst_mod)       
            dest_lst_myd = gdal.Open(filename_lst_myd)              

            Array_angle_mod = dest_angle_mod.GetRasterBand(1).ReadAsArray()
            Array_angle_myd = dest_angle_myd.GetRasterBand(1).ReadAsArray()        
            Array_time_mod = dest_time_mod.GetRasterBand(1).ReadAsArray()        
            Array_time_myd = dest_time_myd.GetRasterBand(1).ReadAsArray()        
            Array_lst_mod = dest_lst_mod.GetRasterBand(1).ReadAsArray()
            Array_lst_myd = dest_lst_myd.GetRasterBand(1).ReadAsArray()        
            
            LST = Array_lst_mod
            Time = Array_time_mod
            LST = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_lst_myd, LST)
            Time = np.where(np.abs(Array_angle_myd)<np.abs(Array_angle_mod), Array_time_myd, Time)        
            
            proj_ex = dest_angle_mod.GetProjection()
            geo_ex = dest_angle_mod.GetGeoTransform()       
            
            DC.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)
            DC.Save_as_tiff(Time_file, Time, geo_ex, proj_ex) 
            
    return()
        
def Calc_Composite_METEO(input_format, Date, example_file, method = "mean", DEM_file = None, lapse = -0.006):
    
    dates_comp = pd.date_range(Date, Date + pd.DateOffset(days=7), freq = "D")
    i = 0
    
    for date_comp in dates_comp:
        print("Find %s" %input_format.format(yyyy=date_comp.year, mm=date_comp.month, dd=date_comp.day))
        filename = glob.glob(input_format.format(yyyy=date_comp.year, mm=date_comp.month, dd=date_comp.day))[0]  
        
        if DEM_file != None:        
            Array = lapse_rate_temp(filename, DEM_file, lapse)
            proj_filename = DEM_file
        else:
            dest_one = gdal.Open(filename)
            Array = dest_one.GetRasterBand(1).ReadAsArray()
            Array[Array==-9999] = np.nan
            proj_filename = filename
        
        if date_comp == dates_comp[0]:
            dest_one = gdal.Open(proj_filename)
            Array_end = np.ones([8, dest_one.RasterYSize, dest_one.RasterXSize])
            geo = dest_one.GetGeoTransform()
            proj = dest_one.GetProjection()
        
        Array_end[i, :, :] = Array
        i+=1

    if method == "mean":
        Array_end = np.nanmean(Array_end, axis = 0)
    if method == "max":
        Array_end = np.nanmax(Array_end, axis = 0)        
    if method == "min":
        Array_end = np.nanmin(Array_end, axis = 0)
    if method == "sum":
        Array_end = np.nansum(Array_end, axis = 0)
        
    dest_mem = DC.Save_as_MEM(Array_end, geo, proj)
    dest_rep = RC.reproject_dataset_example(dest_mem, example_file, 6)    
    Array_end = dest_rep.GetRasterBand(1).ReadAsArray()
    
    return(Array_end)
    
    
def Combine_LST_composite(folders_input_RAW, Startdate, Enddate):
    
    Dates = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)
    
    output_folder_end = os.path.join(folders_input_RAW, "MODIS", "LST", "8_Daily")
    
    if not os.path.exists(output_folder_end):
        os.makedirs(output_folder_end)
    
    for Date in Dates:
        
        print("Create LST composite %s" %Date)
        
        LST_file = os.path.join(output_folder_end, "LST_MCD11A2_K_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
        Time_file = os.path.join(output_folder_end, "Time_MCD11A2_hour_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))    

        if not os.path.exists(LST_file):
            format_lst_mod = os.path.join(folders_input_RAW, "LST", "MOD11", "8_Daily", "LST_MOD11A2_K_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            format_lst_myd = os.path.join(folders_input_RAW, "LST", "MYD11", "8_Daily", "LST_MYD11A2_K_8-daily_%d.%02d.%02d.tif" %(Date.year, Date.month, Date.day))
            os.chdir(os.path.join(folders_input_RAW, "LST", "MOD11", "8_Daily"))
            filename_lst_mod = glob.glob(format_lst_mod)[0]
            os.chdir(os.path.join(folders_input_RAW, "LST", "MYD11", "8_Daily"))           
            filename_lst_myd = glob.glob(format_lst_myd)[0]
                
            dest_lst_mod = gdal.Open(filename_lst_mod)       
            dest_lst_myd = RC.reproject_dataset_example(filename_lst_myd, dest_lst_mod)        
     
            Array_lst_mod = dest_lst_mod.GetRasterBand(1).ReadAsArray()
            Array_lst_myd = dest_lst_myd.GetRasterBand(1).ReadAsArray()        
            
            LST = Array_lst_mod
            Time = np.where(np.isnan(LST), 13.5, 10.5)
            if Array_lst_mod.shape == Array_lst_myd.shape:
                LST = np.where(np.isnan(LST), Array_lst_myd, LST)
                Time = np.where(np.isnan(LST), np.nan, Time)       
            else:
                print("Arrays are not same size for date: %s" %Date)
                
            proj_ex = dest_lst_mod.GetProjection()
            geo_ex = dest_lst_mod.GetGeoTransform()       
        
            DC.Save_as_tiff(LST_file, LST, geo_ex, proj_ex)
            DC.Save_as_tiff(Time_file, Time, geo_ex, proj_ex) 
        else:
            print("LST composite %s already exists" %Date)
            
    return()     

def Thermal_Sharpening_Linear_Forced(surface_temp_up, NDVI_up, NDVI, Box, dest_up, ndvi_fileName, dest_down, watermask = False):

    
    if watermask != False:
       NDVI_up[watermask==1]=np.nan
       surface_temp_up[watermask==1]=np.nan
   
    Buffer_area = int((Box-1)/2) 
    NDVI_up_box =np.empty((Box**3, len(NDVI_up),len(NDVI_up[1]))) * np.nan
    LST_up_box =np.empty((Box**3, len(NDVI_up),len(NDVI_up[1]))) * np.nan
    NDVI_up_box[0, :,:] = NDVI_up
    LST_up_box[0, :,:] = surface_temp_up  
   
    i = 1
    for ypixel in range(0,Buffer_area + 1):

        for xpixel in range(1,Buffer_area + 1):

           if ypixel==0:
                for xpixel in range(1,Buffer_area + 1):
                    NDVI_up_box[int(i),:,0:-xpixel] = NDVI_up[:,xpixel:]
                    LST_up_box[int(i),:,0:-xpixel] = surface_temp_up[:,xpixel:]
                    
                    NDVI_up_box[int(i+1),:,xpixel:] = NDVI_up[:,:-xpixel]
                    LST_up_box[int(i+1),:,xpixel:] = surface_temp_up[:,:-xpixel]
                    i += 2
                    
                for ypixel in range(1,Buffer_area + 1):

                    NDVI_up_box[int(i), ypixel:,:] = NDVI_up[:-ypixel,:]
                    LST_up_box[int(i), ypixel:,:] = surface_temp_up[:-ypixel,:]
                    
                    NDVI_up_box[int(i+1),0:-ypixel,:] = NDVI_up[ypixel:,:]
                    LST_up_box[int(i+1),0:-ypixel,:] = surface_temp_up[ypixel:,:]                    

                    i += 2
                    
                    
                    
           else:
               NDVI_up_box[int(i),0:-xpixel,ypixel:] = NDVI_up[xpixel:,:-ypixel]
               NDVI_up_box[int(i+1),xpixel:,ypixel:] = NDVI_up[:-xpixel,:-ypixel]
               NDVI_up_box[int(i+2),0:-xpixel,0:-ypixel] = NDVI_up[xpixel:,ypixel:]
               NDVI_up_box[int(i+3),xpixel:,0:-ypixel] = NDVI_up[:-xpixel,ypixel:]

               LST_up_box[int(i),0:-xpixel,ypixel:] = surface_temp_up[xpixel:,:-ypixel]
               LST_up_box[int(i+1),xpixel:,ypixel:] = surface_temp_up[:-xpixel,:-ypixel]
               LST_up_box[int(i+2),0:-xpixel,0:-ypixel] = surface_temp_up[xpixel:,ypixel:]
               LST_up_box[int(i+3),xpixel:,0:-ypixel] = surface_temp_up[:-xpixel,ypixel:]
                    
               i += 4

    # Calculate coefficients
    NDVI_up_low = np.nanpercentile(NDVI_up_box, 25, axis = (0))
    NDVI_up_high = np.nanpercentile(NDVI_up_box, 75, axis = (0))
    LST_up_low = np.nanpercentile(LST_up_box, 25, axis = (0))
    LST_up_high = np.nanpercentile(LST_up_box, 75, axis = (0))
    
    # 
    CoefA = (LST_up_low-LST_up_high)/(NDVI_up_high-NDVI_up_low)
    CoefA = CoefA.clip(-30, 0)
    CoefB = LST_up_high - NDVI_up_low * CoefA

    # Define the shape of the surface temperature with the resolution of 400m
    proj = dest_up.GetProjection()
    geo = dest_up.GetGeoTransform()
    
    # Save the coefficients
    CoefA_Downscale = DC.Save_as_MEM(CoefA, geo, proj)
    CoefB_Downscale = DC.Save_as_MEM(CoefB, geo, proj)

    # Downscale the fitted coefficients
    CoefA_Downscale = RC.reproject_dataset_example(CoefA_Downscale, dest_down, 2)    
    CoefB_Downscale = RC.reproject_dataset_example(CoefB_Downscale, dest_down, 2)
    CoefA = CoefA_Downscale.GetRasterBand(1).ReadAsArray()
    CoefB = CoefB_Downscale.GetRasterBand(1).ReadAsArray()

    # Calculate the surface temperature based on the fitted coefficents and NDVI
    temp_surface_sharpened=CoefA*NDVI+CoefB
    temp_surface_sharpened[temp_surface_sharpened < 250] = np.nan
    temp_surface_sharpened[temp_surface_sharpened > 400] = np.nan    

    return(temp_surface_sharpened) 
    



def Thermal_Sharpening_Linear(surface_temp_up, NDVI_up, NDVI, Box, dest_up, ndvi_fileName, dest_down, watermask = False):

    # Creating arrays to store the coefficients
    CoefA=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
    CoefB=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))

    # Fit a second polynominal fit to the NDVI and Thermal data and save the coefficients for each pixel
    # NOW USING FOR LOOPS PROBABLY NOT THE FASTEST METHOD
    for i in range(0,len(surface_temp_up)):
        for j in range(0,len(surface_temp_up[1])):
            if np.isnan(np.sum(surface_temp_up[i,j]))==False and np.isnan(np.sum(NDVI_up[i,j]))==False:
                x_data = NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])))]
                y_data = surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])))]
                if not watermask is False:
                    wm_data = watermask[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])))]
                    x_data = x_data[wm_data==0]
                    y_data = y_data[wm_data==0]
                x_data[~np.isnan(x_data)]
                y_data[~np.isnan(y_data)]
                if len(x_data)>1:
                    coefs = poly.polyfit(x_data, y_data, 1)
                    CoefA[i,j] = coefs[1]
                    CoefB[i,j] = coefs[0]
 
                else:
                    CoefA[i,j] = np.nan
                    CoefB[i,j] = np.nan
            else:
                CoefA[i,j] = np.nan
                CoefB[i,j] = np.nan

    # Define the shape of the surface temperature with the resolution of 400m
    proj = dest_up.GetProjection()
    geo = dest_up.GetGeoTransform()
    
    # Save the coefficients
    CoefA_Downscale = DC.Save_as_MEM(CoefA, geo, proj)
    CoefB_Downscale = DC.Save_as_MEM(CoefB, geo, proj)

    # Downscale the fitted coefficients
    CoefA_Downscale = RC.reproject_dataset_example(CoefA_Downscale, dest_down, 2)    
    CoefB_Downscale = RC.reproject_dataset_example(CoefB_Downscale, dest_down, 2)
    CoefA = CoefA_Downscale.GetRasterBand(1).ReadAsArray()
    CoefB = CoefB_Downscale.GetRasterBand(1).ReadAsArray()

    # Calculate the surface temperature based on the fitted coefficents and NDVI
    temp_surface_sharpened=CoefA*NDVI+CoefB
    temp_surface_sharpened[temp_surface_sharpened < 250] = np.nan
    temp_surface_sharpened[temp_surface_sharpened > 400] = np.nan

    return(temp_surface_sharpened) 

def Thermal_Sharpening(surface_temp_up, NDVI_up, NDVI, Box, dest_up, output_folder, ndvi_fileName, dest_down, watermask = False):

    # Creating arrays to store the coefficients
    CoefA=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
    CoefB=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))
    CoefC=np.zeros((len(surface_temp_up),len(surface_temp_up[1])))

    # Fit a second polynominal fit to the NDVI and Thermal data and save the coefficients for each pixel
    # NOW USING FOR LOOPS PROBABLY NOT THE FASTEST METHOD
    for i in range(0,len(surface_temp_up)):
        for j in range(0,len(surface_temp_up[1])):
            if np.isnan(np.sum(surface_temp_up[i,j]))==False and np.isnan(np.sum(NDVI_up[i,j]))==False:
                x_data = NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])))]
                y_data = surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])))]
                if not watermask is False:
                    wm_data = watermask[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)), int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))][np.logical_and(np.logical_not(np.isnan(NDVI_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]),j + (Box - 1) / 2 + 1))])), np.logical_not(np.isnan(surface_temp_up[int(np.maximum(0, i - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up), i + (Box - 1) / 2 + 1)),int(np.maximum(0, j - (Box - 1) / 2)):int(np.minimum(len(surface_temp_up[1]), j + (Box - 1) / 2 + 1))])))]
                    x_data = x_data[wm_data==0]
                    y_data = y_data[wm_data==0]
                x_data[~np.isnan(x_data)]
                y_data[~np.isnan(y_data)]
                if len(x_data)>6:
                    coefs = poly.polyfit(x_data, y_data, 2)
                    CoefA[i,j] = coefs[2]
                    CoefB[i,j] = coefs[1]
                    CoefC[i,j] = coefs[0]
                else:
                    CoefA[i,j] = np.nan
                    CoefB[i,j] = np.nan
                    CoefC[i,j] = np.nan
            else:
                CoefA[i,j] = np.nan
                CoefB[i,j] = np.nan
                CoefC[i,j] = np.nan

    # Define the shape of the surface temperature with the resolution of 400m
    proj = dest_up.GetProjection()
    geo = dest_up.GetGeoTransform()
    
    # Save the coefficients
    CoefA_Downscale = DC.Save_as_MEM(CoefA, geo, proj)
    CoefB_Downscale = DC.Save_as_MEM(CoefB, geo, proj)
    CoefC_Downscale = DC.Save_as_MEM(CoefC, geo, proj)

    # Downscale the fitted coefficients
    CoefA_Downscaled = RC.reproject_dataset_example(CoefA_Downscale, dest_down)            
    CoefB_Downscaled = RC.reproject_dataset_example(CoefB_Downscale, dest_down)                       
    CoefC_Downscaled = RC.reproject_dataset_example(CoefC_Downscale, dest_down) 
                                        
    CoefA = CoefA_Downscaled.GetRasterBand(1).ReadAsArray()
    CoefB = CoefB_Downscaled.GetRasterBand(1).ReadAsArray()
    CoefC = CoefC_Downscaled.GetRasterBand(1).ReadAsArray()

    # Calculate the surface temperature based on the fitted coefficents and NDVI
    temp_surface_sharpened=CoefA * NDVI**2 + CoefB * NDVI + CoefC
    temp_surface_sharpened[temp_surface_sharpened < 250] = np.nan
    temp_surface_sharpened[temp_surface_sharpened > 400] = np.nan

    return(temp_surface_sharpened)
    
def download_file_from_google_drive(id, destination):
    
    URL = "https://docs.google.com/uc?export=download"

    session = requests.Session()

    response = session.get(URL, params = { 'id' : id }, stream = True)
    token = get_confirm_token(response)

    if token:
        params = { 'id' : id, 'confirm' : token }
        response = session.get(URL, params = params, stream = True)

    save_response_content(response, destination)            

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)          
 