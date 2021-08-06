# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:22:19 2021

@author: timhe
"""

import os
import datetime
from osgeo import gdal
import glob
import pandas as pd
import numpy as np

import watertools
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC

import pywapor

def main(output_folder, Startdate, Enddate, latlim, lonlim, MODIS_WAPOR_folder_in, MODIS_WAPOR_folder_out, LandCover = "GlobCover", Short_Downwards_Radiation = "MSGCCP", RAW_folder = None):

    ############################ Get General inputs ###############################
    
    # Define the input folders
    folders_input_RAW = os.path.join(output_folder, "RAW_composite")
    folder_input_ETLook = os.path.join(output_folder, "ETLook_input_composite_MODIS")
  
    if RAW_folder != None:
        folders_input_RAW = RAW_folder       
  
    # Create folders if not exists
    if not os.path.exists(folders_input_RAW):
        os.makedirs(folders_input_RAW)
    if not os.path.exists(folder_input_ETLook):
        os.makedirs(folder_input_ETLook)
       
    # Define the dates    #!!! is now 8 daily, needs to be changed into dekadal
    Dates = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)
    Timestep_composite = 8
        
    ########################### Download CHIRPS data ############################## 

    ######################## Download Rainfall Data ###############################
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
        pywapor.pre_et_look.download_file_from_google_drive("1pqZnCn-1xkUC7o1csG24hwg22fV57gCH", T_amplitude_global_temp_filename)

    ############### Loop over days for the dynamic data ###############################
    
    # Create the inputs of MODIS
    for Date in Dates:
        
        # try:
        # Define output folder
        folder_input_ETLook_Date = os.path.join(folder_input_ETLook, "%d%02d%02d" %(Date.year, Date.month, Date.day))
        if not os.path.exists(folder_input_ETLook_Date):
            os.makedirs(folder_input_ETLook_Date)

        folder_input_ETLook_Static = os.path.join(folder_input_ETLook, "Static")
        if not os.path.exists(folder_input_ETLook_Static):
            os.makedirs(folder_input_ETLook_Static)


        ##########################################################################
        # Create Albedo, NDVI and SE root
        ##########################################################################    
        
        NDVI_file = os.path.join(folder_input_ETLook_Date, "NDVI_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))            
        if not os.path.exists(NDVI_file):
            
            input_format_NDVI = os.path.join(MODIS_WAPOR_folder_in, "{yyyy}{mm:02d}{dd:02d}", "NDVI_{yyyy}{mm:02d}{dd:02d}.tif")
            NDVI_comp, geo, proj = Calc_Composite_Indicators(input_format_NDVI, Date, Timestep_composite, method = "mean")
            DC.Save_as_tiff(NDVI_file, NDVI_comp, geo, proj)
        
        ALBEDO_file = os.path.join(folder_input_ETLook_Date, "ALBEDO_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
        if not os.path.exists(ALBEDO_file):
            
            input_format_ALBEDO = os.path.join(MODIS_WAPOR_folder_in, "{yyyy}{mm:02d}{dd:02d}", "ALBEDO_{yyyy}{mm:02d}{dd:02d}.tif")
            ALBEDO_comp, geo, proj = Calc_Composite_Indicators(input_format_ALBEDO, Date, Timestep_composite, method = "mean")
            DC.Save_as_tiff(ALBEDO_file, ALBEDO_comp, geo, proj)
                    
        SE_root_file = os.path.join(folder_input_ETLook_Date, "se_root_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))
        if not os.path.exists(SE_root_file):
            
            input_format_SE_root = os.path.join(MODIS_WAPOR_folder_out, "{yyyy}{mm:02d}{dd:02d}", "se_root_{yyyy}{mm:02d}{dd:02d}.tif")
            SE_root_comp, geo, proj = Calc_Composite_Indicators(input_format_SE_root, Date, Timestep_composite, method = "mean")
            DC.Save_as_tiff(SE_root_file, SE_root_comp, geo, proj)

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

        ########################## Create Time rasters ################################

        # Get DOY
        date_middle = Date + datetime.timedelta(days=int(Timestep_composite/2))             
        doy = int(date_middle.strftime("%j"))

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

            input_format = os.path.join(folder_RAW_file_P, "P_CHIRPS.v2.0_mm-day-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif")
            print("LALA", input_format)
            P = pywapor.pre_et_look.Calc_Composite_METEO([input_format], Date, NDVI_file, "mean")
            P[np.isnan(P)] = 0
            DC.Save_as_tiff(P_file, P, geo_ex, proj_ex)                   


        ############################# Download METEO ##################################
        
        # Define the startdates for the METEO
        Date_end = Date + datetime.timedelta(days = int(Timestep_composite - 1)) 
        StartTime = datetime.datetime(Date.year, Date.month, Date.day, 0, 0)
        EndTime = datetime.datetime(Date_end.year, Date_end.month, Date_end.day, 23, 59)                                  
        
        # Download METEO data
        #if Date < datetime.datetime(2016,1,1): 
            # Test https://opendap.nccs.nasa.gov/dods/GEOS-5/MERRAero/hourly/tavg3hr_2d_asm_Nx.ascii?u10m[20:1:24][227:1:255][384:1:420]
        if Date < datetime.datetime(2000,1,1): 
            watertools.Collect.MERRA.daily(folders_input_RAW, ['u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim)
            watertools.Collect.MERRA.daily(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
            str_METEO = "MERRA"
            day_name = "daily"
                
        #elif (Date >= datetime.datetime(2016,1,1) and Date < datetime.datetime(2017,12,1)):     
        elif (Date >= datetime.datetime(2000,1,1) and Date < datetime.datetime(2017,12,1)):     

            watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp'],StartTime, EndTime, latlim, lonlim)
            watertools.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
            str_METEO = "MERRA"   
            day_name = "daily_MERRA2"
            
        else:
            watertools.Collect.GEOS.daily(folders_input_RAW, ['u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp'], StartTime, EndTime, latlim, lonlim)
            watertools.Collect.GEOS.daily(folders_input_RAW, ['t2m'], StartTime, EndTime, latlim, lonlim, data_type = ["mean", "min", "max"])
            str_METEO = "GEOS"
            day_name = "daily"
            
        # Air pressure
        pair_24_0_file = os.path.join(folder_input_ETLook_Date, "Pair_24_0_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))        

        if not os.path.exists(pair_24_0_file):   
            
            folder_RAW_file_pair_24_0 = os.path.join(folders_input_RAW, str_METEO, "Sea_Level_Pressure", day_name) 
            
            input_format_pair_24_sea = os.path.join(folder_RAW_file_pair_24_0, "slp_%s_kpa_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))   

            Pair_24_sea = pywapor.pre_et_look.Calc_Composite_METEO([input_format_pair_24_sea], Date, NDVI_file, "mean")
            Pair_24_sea[np.isnan(Pair_24_sea)] = 0     
            
            DC.Save_as_tiff(pair_24_0_file, Pair_24_sea, geo_ex, proj_ex)
            
        # Specific Humidity
        qv_24_file = os.path.join(folder_input_ETLook_Date, "qv_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
                    
        if not os.path.exists(qv_24_file):   
            folder_RAW_file_qv_24 = os.path.join(folders_input_RAW, str_METEO, "Specific_Humidity", day_name) 
            if str_METEO == "MERRA":
                para = "q2m"
            else:
                para = "qv2m"

            input_format_qv_24 = os.path.join(folder_RAW_file_qv_24, "%s_%s_kg-kg-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(para, str_METEO))

            qv_24 = pywapor.pre_et_look.Calc_Composite_METEO([input_format_qv_24], Date, NDVI_file, "mean")
            qv_24[np.isnan(qv_24)] = 0                       
            
            DC.Save_as_tiff(qv_24_file, qv_24, geo_ex, proj_ex)

        # Air temperature
        Tair_24_file = os.path.join(folder_input_ETLook_Date, "tair_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
        Tair_max_24_file = os.path.join(folder_input_ETLook_Date, "tair_max_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day)) 
        Tair_min_24_file = os.path.join(folder_input_ETLook_Date, "tair_min_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day)) 
                    
        if not (os.path.exists(Tair_24_file) and os.path.exists(Tair_max_24_file) and os.path.exists(Tair_min_24_file)):
            folder_RAW_file_tair_24 = os.path.join(folders_input_RAW, str_METEO, "Air_Temperature", day_name) 
            
            input_format_tair_24 = os.path.join(folder_RAW_file_tair_24, "t2m_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
            input_format_tair_min_24 = os.path.join(folder_RAW_file_tair_24,  "min", "t2mmin_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
            input_format_tair_max_24 = os.path.join(folder_RAW_file_tair_24,  "max", "t2mmax_%s_K_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))        

            tair_24 = pywapor.pre_et_look.Calc_Composite_METEO([input_format_tair_24], Date, NDVI_file, "max", DEM_file = DEM_file, lapse = -0.006)
            tair_24[np.isnan(tair_24)] = 0                       
            if np.nanmax(tair_24>270):
                tair_24 = tair_24 -273.15
                
            DC.Save_as_tiff(Tair_24_file, tair_24, geo_ex, proj_ex)

            tair_min_24 = pywapor.pre_et_look.Calc_Composite_METEO([input_format_tair_min_24], Date, NDVI_file, "max", DEM_file = DEM_file, lapse = 0.0)
            tair_min_24[np.isnan(tair_min_24)] = 0                       
            if np.nanmax(tair_min_24>270):
                tair_min_24 = tair_min_24 -273.15
                                    
            DC.Save_as_tiff(Tair_min_24_file, tair_min_24, geo_ex, proj_ex)                   
                
            tair_max_24 = pywapor.pre_et_look.Calc_Composite_METEO([input_format_tair_max_24], Date, NDVI_file, "max", DEM_file = DEM_file, lapse = 0.0)
            tair_max_24[np.isnan(tair_max_24)] = 0                       
            if np.nanmax(tair_max_24>270):
                tair_max_24 = tair_max_24 -273.15                       
            DC.Save_as_tiff(Tair_max_24_file, tair_max_24, geo_ex, proj_ex)                   
    
        # Wind Speed
        wind_24_file = os.path.join(folder_input_ETLook_Date, "wind_24_%d%02d%02d.tif" %(Date.year, Date.month, Date.day))  
                    
        if not os.path.exists(wind_24_file):
                    
            folder_RAW_file_u2_24 = os.path.join(folders_input_RAW, str_METEO, "Eastward_Wind", day_name)      
            folder_RAW_file_v2_24 = os.path.join(folders_input_RAW, str_METEO, "Northward_Wind", day_name) 

            input_format_u2_24 = os.path.join(folder_RAW_file_u2_24, "u2m_%s_m-s-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))
            input_format_v2_24 = os.path.join(folder_RAW_file_v2_24, "v2m_%s_m-s-1_daily_{yyyy}.{mm:02d}.{dd:02d}.tif" %(str_METEO))

            u2_24= pywapor.pre_et_look.Calc_Composite_METEO([input_format_u2_24], Date, NDVI_file, "mean")
            u2_24[np.isnan(u2_24)] = 0                       

            v2_24= pywapor.pre_et_look.Calc_Composite_METEO([input_format_v2_24], Date, NDVI_file, "mean")
            v2_24[np.isnan(v2_24)] = 0    
            
            wind_24 = np.sqrt(u2_24**2 + v2_24 **2) 
                
            DC.Save_as_tiff(wind_24_file, wind_24, geo_ex, proj_ex)                      

                
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
                import pywapor.general.landcover_converter as LCC
                
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
                desttairamp = RC.reproject_dataset_example(T_amplitude_global_temp_filename, NDVI_file, method=2)
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
                
            Dates_End_Rads = Date + datetime.timedelta(days = int(Timestep_composite - 1)) 
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
                    destswgnet = RC.reproject_dataset_example(dest_swgnet_mean, NDVI_file, method=2)
                    del geo_trans
                    
                else:
                    destswgnet = RC.reproject_dataset_example(os.path.join(folder_RAW_file_trans, filename_trans), NDVI_file, method=2)
                        
                swgnet_one = destswgnet.GetRasterBand(1).ReadAsArray()
                if i == 0:
                    swgnet = np.ones([8, destswgnet.RasterYSize, destswgnet.RasterXSize])
                swgnet[i, :, :] = swgnet_one
                i+=1
                
            trans = np.nanmean(swgnet,axis = 0) / Ra24_flat
            DC.Save_as_tiff(Trans_file, trans, geo_ex, proj_ex)


        # except Exception as e:
        #     print("No ETLook input dataset for %s" %Date)
        #     print(e) 
            
    return()


def Calc_Composite_Indicators(input_format, Date, Timestep_composite, method = "mean"):
    
    dates_comp = pd.date_range(Date, Date + datetime.timedelta(days=int(Timestep_composite-1)), freq = "D")
    i = 0
    
    for date_comp in dates_comp:
        print("Find %s" %input_format.format(yyyy=date_comp.year, mm=date_comp.month, dd=date_comp.day))
        filename = glob.glob(input_format.format(yyyy=date_comp.year, mm=date_comp.month, dd=date_comp.day))[0]  
        
        dest_one = gdal.Open(filename)
        Array = dest_one.GetRasterBand(1).ReadAsArray()
        Array[Array==-9999] = np.nan
        proj_filename = filename
    
        if date_comp == dates_comp[0]:
            dest_one = gdal.Open(proj_filename)
            geo = dest_one.GetGeoTransform()
            proj = dest_one.GetProjection()
            Array_end = np.ones([int(Timestep_composite), dest_one.RasterYSize, dest_one.RasterXSize]) * np.nan
        
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
        
    return(Array_end, geo, proj)