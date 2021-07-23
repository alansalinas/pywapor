# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:43:42 2020

@author: timhe
"""

#%%
import re
import os
from osgeo import osr
import glob
from osgeo import gdal
import datetime 
import numpy as np
import pandas as pd
from pyproj import Proj, transform
import time
import shutil
import SEBAL.pySEBAL.pySEBAL_code as sebal
import watertools
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC
import warnings

def Find_nearest_time(Dates, Date):
    return min(Dates, key=lambda x: abs(x - Date))               
            
def Process_LS_Image(LS_file, FOLDERS, filename_DEM, UTM_Zone, Startdate_IrriEngine, GMT_offset):
    
    print("Process: %s" %LS_file) 
    
    # Define Name of Landsat Image without band number
    Name_Landsat_Image = os.path.splitext(LS_file)[0][:-4]
    Path_row = os.path.splitext(LS_file)[0].split("_")[2]
  
    # the path to the MTL file of landsat
    Landsat_meta_fileName = os.path.join(FOLDERS[0], '%s_MTL.txt' % Name_Landsat_Image)

    # read out the general info out of the MTL file
    print("collect time information from the MTL file: %s" %Landsat_meta_fileName) 
    year, DOY, hour_GTM, minutes_GTM, UTM_Zone_LS, Sun_elevation = info_general_metadata(Landsat_meta_fileName)    

    # Get startdate IrriEngine
    Date_Start = datetime.datetime.strptime(Startdate_IrriEngine, "%Y-%m-%d") + datetime.timedelta(days=1)
    Date_End = datetime.datetime.strptime(Startdate_IrriEngine, "%Y-%m-%d") - pd.DateOffset(months=3)
    Date_Start_or = Date_Start.toordinal()
    Date_End_or = Date_End.toordinal()
    
    # Get date LS in ordinal
    Date_LS = datetime.datetime.strptime("%d%03d" %(year, DOY), "%Y%j")
    Date_LS_or = Date_LS.toordinal()    
    
    if (Date_LS_or < Date_Start_or and Date_LS_or > Date_End_or):

        # Define Time parameter using year and DOY
        TIME = datetime.datetime.strptime("%d%d_%02d%02d"%(year, DOY, hour_GTM, minutes_GTM), "%Y%j_%H%M")    
        
        # Define output cloud
        filename_Cloud = os.path.join(FOLDERS[7], "LS_CLOUD_%s_%d%02d%02d.tif" %(Path_row, TIME.year, TIME.month, TIME.day))  
     
        # Check if cloud mask already exists
        if not os.path.exists(filename_Cloud):
                   
            # Get landsat number
            Landsat_nr = int(LS_file.split("_")[0][2:4])
            
            if 'shape_dem' not in locals():
                            
                print("open dem file: %s" %filename_DEM) 
                dest_dem, DEM_resh, geo, proj, ncol, nrow, shape_dem = Create_Example_UTM_File(filename_DEM, 30, UTM_Zone)
                UTM_Zone = RC.Get_epsg(dest_dem)            
                
                # Create Latitude and Longitude files in UTM
                lat_proy_UTM, lon_proy_UTM  = DEM_lat_lon(dest_dem)      
            
                # Calculate Latitude and Longitude files in WGS84
                input_projection = Proj(init="epsg:%s" %UTM_Zone)
                output_projection = Proj(init="epsg:4326")
                lon_proy, lat_proy = transform(input_projection, output_projection, lon_proy_UTM.flatten(), lat_proy_UTM.flatten())
                lon_proy.resize((nrow, ncol))
                lat_proy.resize((nrow, ncol))
     
                # Calculate aspect and slope
                slope, aspect = Calc_Gradient(DEM_resh, 30)[2:]
    
            # calculate GTM hour
            hour_loc = hour_GTM + GMT_offset
            minutes_loc = minutes_GTM
    
            # Calculate Cosinus zenith angle
            print("calculate the cosinus zenith angle") 
            cos_zn, dr = Calc_Coz_zn_dr(GMT_offset, DOY, hour_loc, minutes_loc, lon_proy, lat_proy, slope, aspect)
    
            # Define bands used for each Landsat number
            if Landsat_nr == 5 or Landsat_nr == 7:
                Bands = np.array([1, 2, 3, 4, 5, 7, 6])
            elif Landsat_nr == 8:
               Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
            else:
                raise RuntimeError('Landsat image not supported, use Landsat 7 or 8')           
    
            # Get correction parameters LS
            
            # Open MTL landsat and get the correction parameters
            print("collect band specific information") 
            Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)
    
            # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
            # for the different Landsat images (L5, L7, or L8)
            ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
            ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
            ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])   
            
            # Create MASK for no data values LS
            filename_BQA = os.path.join(FOLDERS[0], '%s_BQA.TIF' % (Name_Landsat_Image))
            LS_BQA = Open_landsat(filename_BQA, dest_dem)
            LS_BQA = np.where(LS_BQA==0, 1, LS_BQA)
            
            # Cloud Thresholds LS
            if Landsat_nr == 8:
                Cloud_Treshold = np.array([1,2,2722,2720,2724,2728,2732,3744,3748,3752,3756])
                Real_Clouds = np.array([2800,2804,2808,2812,6816,6820,6824,6828,6848,6852,6856,6860,6896,6900,6904,6908,7072,7076,7080,7084,7104,7108,7112,7116,7840,7844,7848,7852,7872,7872,7876,7880,7884]) 
            if Landsat_nr == 5 or Landsat_nr == 7:
                Cloud_Treshold = np.array([1,672,676,680,684,1696,1700,1704,1708])
                Real_Clouds = np.array([752,756,760,764]) 
               
            if Landsat_nr == 8:                
                Cloud_Buffer_Area = np.where(np.isin(LS_BQA, Real_Clouds), 0, 1)
            if Landsat_nr == 5 or Landsat_nr == 7:
                Cloud_Buffer_Area = np.where(np.isin(LS_BQA, Real_Clouds), 0, 1)
            Cloud_Buffer_Area = RC.Create_Buffer(Cloud_Buffer_Area, 3)
            Cloud_Buffer_Area = np.where(Cloud_Buffer_Area==0, 1, 0)
                
            print("create cloud array")     
            QC_mask_Cloud_first = np.where(np.isin(LS_BQA, Cloud_Treshold), 1, np.nan)
            
            # Apply buffer around pixels
            print("apply buffer of 15 pixels of 30m around clouds")                 
            QC_mask_Cloud_second = RC.Create_Buffer(Cloud_Buffer_Area, 15)
            QC_mask_Cloud_second = np.where(QC_mask_Cloud_second == 1, np.nan, 1)

            QC_mask_Cloud = (QC_mask_Cloud_first + QC_mask_Cloud_second)/2
            QC_mask_Cloud[LS_BQA==1] = np.nan
                       
            # Create 3D array to store Spectral radiance and Reflectivity for each band
            print("get reflectance and spectral radiance of the bands")              
            Reflect, Spec_Rad = Landsat_Reflect(Bands, FOLDERS[0], Name_Landsat_Image, shape_dem, QC_mask_Cloud, Lmax, Lmin, ESUN_L5, ESUN_L7, ESUN_L8, cos_zn, dr, Landsat_nr, dest_dem)
            
            # Calculate NDVI
            print("calculate spectral indices")                  
            NDVI = Calc_NDVI(Reflect)
            NDII = (Reflect[:, :, 3] - Reflect[:, :, 4])/(Reflect[:, :, 3] + Reflect[:, :, 4])
            VSDI = 1-((Reflect[:, :, 4] - Reflect[:, :, 0]) + (Reflect[:, :, 2] - Reflect[:, :, 0]))
        
            # Surface albedo:
            Surf_albedo = (0.3 * Reflect[:, :, 0] + 0.277 * Reflect[:, :, 1] +
                           0.233 * Reflect[:, :, 2] + 0.143 * Reflect[:, :, 3] +
                           0.036 * Reflect[:, :, 4] + 0.012 * Reflect[:, :, 5] -
                           0.03 ) / np.power(0.89 , 2)            
            
            # B12 of S2
            B12 = Reflect[:, :, 5]
            
            # Remove bad pixels
            NDVI[np.isnan(QC_mask_Cloud)] = np.nan
            NDII[np.isnan(QC_mask_Cloud)] = np.nan
            VSDI[np.isnan(QC_mask_Cloud)] = np.nan
            Surf_albedo[np.isnan(QC_mask_Cloud)] = np.nan
            B12[np.isnan(QC_mask_Cloud)] = np.nan    
            
            # Guess the first watermask
            water_mask_temp = np.where(NDVI < 0.0, 1.0, 0)
            
            # Calculate LAI
            FPAR, tir_emis, Nitrogen, vegt_cover, LAI, b10_emissivity = Calc_vegt_para(NDVI, water_mask_temp)
    
            # Get meteo data
            print("get meteo data from excel file")        

            RH_format_inst = os.path.join(FOLDERS[1], "Weather_Data", "Model", "GLDAS", "three_hourly", "rh_f_inst", "Hum_GLDAS-NOAH_percentage_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
            temp_format_inst = os.path.join(FOLDERS[1], "Weather_Data", "Model", "GLDAS", "three_hourly", "tair_f_inst", "Tair_GLDAS-NOAH_C_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
             
            dest_rh_inst = RC.reproject_dataset_example(RH_format_inst.format(yyyy=TIME.year, mm=TIME.month, dd=TIME.day, HH=hour_take), dest_dem, 2)
            Humidity_inst = dest_rh_inst.GetRasterBand(1).ReadAsArray()
            Humidity_inst[Humidity_inst==0] = np.nan
            Humidity_inst[np.isnan(Humidity_inst)] = np.nanmean(Humidity_inst)
            
            dest_temp_inst = RC.reproject_dataset_example(temp_format_inst.format(yyyy=TIME.year, mm=TIME.month, dd=TIME.day, HH=hour_take), dest_dem, 2)
            Temperature_inst = dest_temp_inst.GetRasterBand(1).ReadAsArray()
            Temperature_inst[Temperature_inst==0] = np.nan
            Temperature_inst[np.isnan(Temperature_inst)] = np.nanmean(Temperature_inst)
           
                      
            try: 
                
                # Calculate vapour pressure
                esat_inst = 0.6108 * np.exp(17.27 * Temperature_inst / (Temperature_inst + 237.3)) 
                eact_inst = Humidity_inst * esat_inst / 100
                
                # Define amount of bands
                if Landsat_nr == 8:
                    Bands_thermal = 2
                else:
                    Bands_thermal = 1
                
                # Calculate Surface temperature
                print("calculate lst of landsat")               
                therm_data = Landsat_therm_data(Bands, FOLDERS[0], Name_Landsat_Image, shape_dem, QC_mask_Cloud, dest_dem)      
                therm_data[therm_data==0] = np.nan
                Surface_temp, cloud_mask_temp = Calc_surface_water_temp(Temperature_inst, Landsat_nr, Lmax, Lmin, therm_data, b10_emissivity, k1_c, k2_c, eact_inst, shape_dem, water_mask_temp, Bands_thermal, 0.91 , 0.866 , 30, 1)
          
                # Define output names
                LS_fileT = os.path.join(FOLDERS[8], "LS_LST_%s_%d%02d%02d_%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day, TIME.hour, TIME.minute))
                DC.Save_as_tiff(LS_fileT, Surface_temp, geo, proj)   
                
            except:
                print("Temperature is negative, so LS %s is not processed" %Name_Landsat_Image)
                
               
            Surf_albedo[Surf_albedo>1] = np.nan
            NDVI[NDVI>1] = np.nan
            VSDI[VSDI==0.] = np.nan
            VSDI[VSDI==1] = np.nan
            NDII[NDII==0.] = np.nan 

            filename_NDVI = os.path.join(FOLDERS[2], "LS_NDVI_%s_%d%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day))
            filename_VSDI = os.path.join(FOLDERS[4], "LS_VSDI_%s_%d%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day))
            filename_NDII = os.path.join(FOLDERS[5], "LS_NDII_%s_%d%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day))
            filename_Albedo = os.path.join(FOLDERS[3], "LS_Albedo_%s_%d%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day))
            filename_B12 = os.path.join(FOLDERS[6], "LS_B12_%s_%d%02d%02d.tif" %(Path_row,TIME.year, TIME.month, TIME.day))
            
            # Save NDVI, Albedo, VSDI, NDII, B12, Cloud
            print("save bands of landsat")                
            DC.Save_as_tiff(filename_Cloud, QC_mask_Cloud, geo, proj)   
            DC.Save_as_tiff(filename_NDVI, NDVI, geo, proj)      
            DC.Save_as_tiff(filename_VSDI, VSDI, geo, proj)      
            DC.Save_as_tiff(filename_NDII, NDII, geo, proj)      
            DC.Save_as_tiff(filename_Albedo, Surf_albedo, geo, proj)      
            DC.Save_as_tiff(filename_B12, B12, geo, proj)   
              
    return()  

#------------------------------------------------------------------------------
def info_general_metadata(filename):
    """
    This function retrieves general information of the Landsat image
    (date and time aquired, UTM zone, sun elevation) from the
    metadata file.

    """
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)SCENE_CENTER_TIME(.*)", line): # search in metadata for line SCENE_CENTER_TIME
            words = line.split()# make groups of words which are divided by an open space
            time_list = words[2].split(':', 2) # Take the second word of words and split the word which are divided by :
            if len(time_list[0])== 3:
                time_list[0]=time_list[0][1:3]
                time_list[2]=time_list[2][0:-1]
            hour = float(time_list[0]) # take the first word of time_list
            minutes = float(time_list[1]) + float(time_list[2][:-1]) / 60 # Take the second and third word of time_list and place :-1 to remove Z behined minutes
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)DATE_ACQUIRED(.*)", line):
            words = line.split()
            DOY = time.strptime(words[2], "%Y-%m-%d").tm_yday
            year = time.strptime(words[2], "%Y-%m-%d").tm_year
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)UTM_ZONE(.*)", line):
            words = line.split()
            UTM_Zone = int(words[2])
    Landsat_meta = open(filename, "r")  # Open metadata file
    for line in Landsat_meta:
        if re.match("(.*)SUN_ELEVATION(.*)", line):
            words = line.split()
            Sun_elevation = float(words[2])

    return year, DOY, hour, minutes, UTM_Zone, Sun_elevation

def Calc_vegt_para(NDVI,water_mask_temp):
    """
    Calculates the Fraction of PAR, Thermal infrared emissivity, Nitrogen, Vegetation Cover, LAI, b10_emissivity
    """
    # Fraction of PAR absorbed by the vegetation canopy (FPAR):
    FPAR = -0.161 + 1.257 * NDVI
    FPAR[NDVI < 0.125] = 0.0

    # Termal infrared emissivity
    tir_emis = 1.009 + 0.047 * np.log(NDVI)
    tir_emis[np.logical_or(water_mask_temp == 1.0, water_mask_temp == 2.0)] = 1.0
    tir_emis[np.logical_and(NDVI < 0.125, water_mask_temp == 0.0)] = 0.92

    # Vegetation Index - Regression model from Bagheri et al. (2013)
    VI = 38.764 * np.square(NDVI) - 24.605 * NDVI + 5.8103

    # Nitrogen computation
    Nitrogen = np.copy(VI)
    Nitrogen[VI <= 0.0] = 0.0
    Nitrogen[NDVI <= 0.0] = 0.0

    # Vegetation cover:
    vegt_cover = 1 - np.power((0.8 - NDVI)/(0.8 - 0.125), 0.7)
    vegt_cover[NDVI < 0.125] = 0.0
    vegt_cover[NDVI > 0.8] = 0.99

    # Leaf Area Index (LAI)
    LAI_1 = np.log(-(vegt_cover - 1)) / -0.45
    LAI_1[LAI_1 > 8] = 8.0
    LAI_2 = (9.519 * np.power(NDVI, 3) + 0.104 * np.power(NDVI, 2) +
             1.236 * NDVI - 0.257)

    LAI = (LAI_1 + LAI_2) / 2.0  # Average LAI
    LAI[LAI < 0.001] = 0.001

    b10_emissivity = np.where(LAI <= 3.0, 0.95 + 0.01 * LAI, 0.98)
    b10_emissivity[water_mask_temp != 0.0] = 1.0

    return(FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity)

def Calc_NDVI(Reflect):
    """
    This function calculates and returns the Surface albedo, NDVI by using the refectance from the landsat image.
    """
    # Computation of Normalized Difference Vegetation Index (NDVI)
    NDVI = ((Reflect[:, :, 3] - Reflect[:, :, 2]) /
            (Reflect[:, :, 3] + Reflect[:, :, 2]))

    return(NDVI)

def Calc_Gradient(dataset,pixel_spacing):
    """
    This function calculates the slope and aspect of a DEM map.
    """
    # constants
    deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
    rad2deg = 180.0 / np.pi  # Factor to transform from rad to degree

    # Calculate slope
    x, y = np.gradient(dataset, pixel_spacing, pixel_spacing)
    hypotenuse_array = np.hypot(x,y)
    slope = np.arctan(hypotenuse_array) * rad2deg
    #slope = np.arctan(np.sqrt(np.square(x/pixel_spacing) + np.square(y/pixel_spacing))) * rad2deg

    # calculate aspect
    aspect = np.arctan2(y/pixel_spacing, -x/pixel_spacing) * rad2deg
    aspect = 180 + aspect

    return(deg2rad,rad2deg,slope,aspect)


def Calc_surface_water_temp(Temp_inst,Landsat_nr,Lmax,Lmin,therm_data,b10_emissivity,k1_c,k2_c,eact,shape_lsc,water_mask_temp,Bands_thermal,Rp,tau_sky,surf_temp_offset,Image_Type):
    """
    Calculates the surface temperature and create a water mask
    """

    # Spectral radiance for termal
    if Landsat_nr == 8:
        if Bands_thermal == 1:
            k1 = k1_c[0]
            k2 = k2_c[0]
            L_lambda_b10 = (Lmax[-1] - Lmin[-1]) / (65535-1) * therm_data[:, :, 0] + Lmin[-1]

            # Get Temperature
            Surface_temp = Get_Thermal(L_lambda_b10,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

        elif Bands_thermal == 2:
            L_lambda_b10 = (Lmax[-2] - Lmin[-2]) / (65535-1) * therm_data[:, :, 0] + Lmin[-2]
            L_lambda_b11 = (Lmax[-1] - Lmin[-1]) / (65535-1) * therm_data[:, :, 1] + Lmin[-1]

            # Brightness temperature
            # From Band 10:
            Temp_TOA_10 = (k2_c[0] / np.log(k1_c[0] / L_lambda_b10 + 1.0))
            # From Band 11:
            Temp_TOA_11 = (k2_c[1] / np.log(k1_c[1] / L_lambda_b11 + 1.0))
            # Combined:
            Surface_temp = (Temp_TOA_10 + 1.378 * (Temp_TOA_10 - Temp_TOA_11) +
                           0.183 * np.power(Temp_TOA_10 - Temp_TOA_11, 2) - 0.268 +
                           (54.30 - 2.238 * eact) * (1 - b10_emissivity))

    elif Landsat_nr == 7:
        k1=666.09
        k2=1282.71
        L_lambda_b6 = (Lmax[-1] - Lmin[-1]) / (256-1) * therm_data[:, :, 0] + Lmin[-1]

        # Brightness temperature - From Band 6:
        Surface_temp = Get_Thermal(L_lambda_b6,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

    elif Landsat_nr == 5:
        k1=607.76
        k2=1260.56
        L_lambda_b6 = ((Lmax[-1] - Lmin[-1]) / (256-1) * therm_data[:, :, 0] +
                       Lmin[-1])

       # Brightness temperature - From Band 6:
        Surface_temp = Get_Thermal(L_lambda_b6,Rp,Temp_inst,tau_sky,b10_emissivity,k1,k2)

    # Surface temperature
    Surface_temp = Surface_temp.clip(230.0, 360.0)

    # Cloud mask:
    temp_water = np.zeros((shape_lsc[1], shape_lsc[0]))
    temp_water = np.copy(Surface_temp)
    temp_water[water_mask_temp == 0.0] = np.nan
    temp_water_sd = np.nanstd(temp_water)     # Standard deviation
    temp_water_mean = np.nanmean(temp_water)  # Mean
    print('Mean water temperature = ', '%0.3f (Kelvin)' % temp_water_mean)
    print('SD water temperature = ', '%0.3f (Kelvin)' % temp_water_sd)
    cloud_mask = np.zeros((shape_lsc[1], shape_lsc[0]))
    #cloud_mask[Surface_temp < np.minimum((temp_water_mean - 1.0 * temp_water_sd -
    #           surf_temp_offset),290)] = 1.0

    return(Surface_temp, cloud_mask)

def Get_Thermal(lambda_b10,Rp,Temp_inst,tau_sky,TIR_Emissivity,k1,k2):

    # Narrow band downward thermal radiation from clear sky, rsky (W/m2/sr/µm)
    rsky = (1.807E-10 * np.power(Temp_inst + 273.15, 4) * (1 - 0.26 *
            np.exp(-7.77E-4 * np.power((-Temp_inst), -2))))
    
    print('Rsky = ', '%0.3f (W/m2/sr/µm)' % np.nanmean(rsky))

    # Corrected thermal radiance from the surface, Wukelikc et al. (1989):
    correc_lambda_b10 = ((lambda_b10 - Rp) / tau_sky -
                               (1.0 - TIR_Emissivity) * rsky)
    
    # Brightness temperature - From Band 10:
    Temp_TOA = (k2 / np.log(TIR_Emissivity * k1 /
                       correc_lambda_b10 + 1.0))

    return(Temp_TOA)

def DEM_lat_lon(data_in):
    """
    This function retrieves information about the latitude and longitude of the
    DEM map.

    """
    if str(type(data_in)) == "<class 'str'>":
        g = gdal.Open(data_in)     # Open DEM
    elif str(type(data_in)) == "<class 'osgeo.gdal.Dataset'>":
        g = data_in
    else:
        raise RuntimeError("Dataformat not suported")
        
    geo_t = g.GetGeoTransform()     # Get the Geotransform vector:
    x_size = g.RasterXSize          # Raster xsize - Columns
    y_size = g.RasterYSize          # Raster ysize - Rows

    # create a longitude and a latitude array
    lon = np.zeros((y_size, x_size))
    lat = np.zeros((y_size, x_size))
    for col in np.arange(x_size):
        lon[:, col] = geo_t[0] + col * geo_t[1] + geo_t[1]/2
        # ULx + col*(E-W pixel spacing) + E-W pixel spacing
    for row in np.arange(y_size):
        lat[row, :] = geo_t[3] + row * geo_t[5] + geo_t[5]/2
        # ULy + row*(N-S pixel spacing) + N-S pixel spacing,
        # negative as we will be counting from the UL corner

    return(lat, lon)

def w_time(GMT,lon_proy, DOY):
    """
    This function computes the hour angle (radians) of an image given the
    local time, longitude, and day of the year.

    """
    nrow, ncol = lon_proy.shape

    # Difference of the local time (LT) from Greenwich Mean Time (GMT) (hours):
    delta_GTM = lon_proy[int(nrow/2), int(ncol/2)] * 24 / 360
    if np.isnan(delta_GTM) == True:
         delta_GTM = np.nanmean(lon_proy) * np.nanmean(lon_proy)  * 24 / 360

    # Local Standard Time Meridian (degrees):
    LSTM = 15 * delta_GTM

    # Ecuation of time (EoT, minutes):
    B = 360./365 * (DOY-81)  # (degrees)
    EoT = 9.87*np.sin(np.deg2rad(2*B))-7.53*np.cos(np.deg2rad(B))-1.5*np.sin(np.deg2rad(B))

    # Net Time Correction Factor (minutes) at the center of the image:
    TC = 4 * (lon_proy - LSTM) + EoT     # Difference in time over the longitude
    LST = GMT + delta_GTM + TC/60         # Local solar time (hours)
    HRA = 15 * (LST-12)                  # Hour angle HRA (degrees)
    w = np.deg2rad(HRA)                    # Hour angle HRA (radians)
    return w

def reproject_dataset(dataset, pixel_spacing, UTM_Zone, fit_extend = False):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection
    """

    # 1) Open the dataset
    try:                               #!!! Doe deze met de type check
        g = gdal.Open(dataset)
    except:
        g = dataset
   
     # Define the EPSG code...
    EPSG_code = int(UTM_Zone)
    epsg_to = int(EPSG_code)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    try:
        proj = g.GetProjection()
        Proj_in=proj.split('EPSG","')
        epsg_from=int((str(Proj_in[-1]).split(']')[0])[0:-1])
    except:
        epsg_from = int(4326)    # Get the Geotransform vector:
    geo_t = g.GetGeoTransform()

    # Vector components:
    # 0- The Upper Left easting coordinate (i.e., horizontal)
    # 1- The E-W pixel spacing
    # 2- The rotation (0 degrees if image is "North Up")
    # 3- The Upper left northing coordinate (i.e., vertical)
    # 4- The rotation (0 degrees)
    # 5- The N-S pixel spacing, negative as it is counted from the UL corner
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize

    epsg_to = int(epsg_to)

    # 2) Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    osng = osr.SpatialReference()
    osng.ImportFromEPSG(epsg_to)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(epsg_from)

    inProj = Proj(init='epsg:%d' %epsg_from)
    outProj = Proj(init='epsg:%d' %epsg_to)

    # Remove a part of image
    nrow_skip = 0 #round((0.06*y_size)/2)
    ncol_skip = 0 #round((0.06*x_size)/2)

    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the to
    ulx, uly = transform(inProj,outProj,geo_t[0] + nrow_skip * geo_t[1], geo_t[3] + nrow_skip * geo_t[5])
    lrx, lry = transform(inProj,outProj,geo_t[0] + geo_t[1] * (x_size-ncol_skip),
                                        geo_t[3] + geo_t[5] * (y_size-nrow_skip))

    # See how using 27700 and WGS84 introduces a z-value!
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')

    if fit_extend == True:       
        ulx = np.ceil(ulx/pixel_spacing) * pixel_spacing + 0.5 * pixel_spacing
        uly = np.floor(uly/pixel_spacing) * pixel_spacing - 0.5 * pixel_spacing
        lrx = np.floor(lrx/pixel_spacing) * pixel_spacing - 0.5 * pixel_spacing
        lry = np.ceil(lry/pixel_spacing) * pixel_spacing + 0.5 * pixel_spacing

    else:       
        ulx = np.ceil(ulx/pixel_spacing) * pixel_spacing
        uly = np.floor(uly/pixel_spacing) * pixel_spacing
        lrx = np.floor(lrx/pixel_spacing) * pixel_spacing
        lry = np.ceil(lry/pixel_spacing) * pixel_spacing
 
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    col = int((lrx - ulx)/pixel_spacing)
    rows = int((uly - lry)/pixel_spacing)

    # Re-define lr coordinates based on whole number or rows and columns
    (lrx, lry) = (ulx + col * pixel_spacing, uly -
                  rows * pixel_spacing)

    dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)

    if dest is None:
        print('input folder to large for memory, clip input map')

   # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2], uly,
               geo_t[4], - pixel_spacing)

    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(osng.ExportToWkt())

    # Perform the projection/resampling
    gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(),gdal.GRA_CubicSpline)

    return dest, ulx, lry, lrx, uly, epsg_to

    
def Create_Example_UTM_File(filename_ex, pixel_spacing, proj):

    # Reproject DEM to fit LS scene
    dest, ulx, lry, lrx, uly, epsg_to = reproject_dataset(
            filename_ex, pixel_spacing, UTM_Zone = proj, fit_extend = True) 
    
    # Find example geotransform info
    Array = dest.GetRasterBand(1).ReadAsArray()
    geo = dest.GetGeoTransform()
    proj = dest.GetProjection()
    ncol = dest.RasterXSize        # Get the reprojected dem column size
    nrow = dest.RasterYSize        # Get the reprojected dem row size
    shape = [ncol, nrow]   
    
    return(dest, Array, geo, proj, ncol, nrow, shape)
    
def Calc_Coz_zn_dr(GMT_offset, DOY, hour_loc, minutes_loc, lon_proy, lat_proy, slope, aspect):
    """
    Calculates the extraterrestiral solar radiation by using the date, slope and aspect.
    """

    # Constants
    deg2rad = np.pi / 180.0  # Factor to transform from degree to rad
    Min_cos_zn = 0.1  # Min value for cos zenith angle
    Max_cos_zn = 1.0  # Max value for cos zenith angle
     
    try:
        GMT_time = float(hour_loc) - GMT_offset + float(minutes_loc)/60  # Local time (hours)
        Loc_time = float(hour_loc) + float(minutes_loc)/60  # Local time (hours)        
    except:
        GMT_time = np.float_(hour_loc) - GMT_offset + np.float_(minutes_loc)/60  # Local time (hours)
        Loc_time = np.float_(hour_loc) + np.float_(minutes_loc)/60  # Local time (hours)

    print('  Local Time: ', '%0.3f' % np.nanmean(Loc_time))
    print('  GMT Time: ', '%0.3f' % np.nanmean(GMT_time))    
    print('  Difference of local time (LT) from Greenwich (GMT): ', GMT_offset)

    # 1. Calculation of extraterrestrial solar radiation for slope and aspect
    # Computation of Hour Angle (HRA = w)
    B = 360./365 * (DOY-81)           # (degrees)
    # Computation of cos(theta), where theta is the solar incidence angle
    # relative to the normal to the land surface
    delta=np.arcsin(np.sin(23.45*deg2rad)*np.sin(np.deg2rad(B))) # Declination angle (radians)
    phi = lat_proy * deg2rad                                     # latitude of the pixel (radians)
    s = slope * deg2rad                                          # Surface slope (radians)
    gamma = (aspect-180) * deg2rad                               # Surface aspect angle (radians)
    w = w_time(GMT_time, lon_proy, DOY)                            # Hour angle (radians)
    a,b,c = Constants(delta,s,gamma,phi)
    cos_zn= AngleSlope(a,b,c,w)
    cos_zn = cos_zn.clip(Min_cos_zn, Max_cos_zn)

    print('Average Cos Zenith Angle: ', '%0.3f (Radians)' % np.nanmean(cos_zn))

    dr = 1 + 0.033 * np.cos(DOY*2*np.pi/365)  # Inverse relative distance Earth-Sun

    return(cos_zn, dr)    

def AngleSlope(a,b,c,w):
    '''
    Based on Richard G. Allen 2006
    Calculate the cos zenith angle by using the hour angle and constants
    '''
    angle = -a + b*np.cos(w) + c*np.sin(w)

    return(angle)

def Constants(delta,s,gamma,phi):
    '''
    Based on Richard G. Allen 2006 equation 11
    determines constants for calculating the exterrestial solar radiation
    '''
    a = np.sin(delta)*np.cos(phi)*np.sin(s)*np.cos(gamma) - np.sin(delta)*np.sin(phi)*np.cos(s)
    b = np.cos(delta)*np.cos(phi)*np.cos(s) + np.cos(delta)*np.sin(phi)*np.sin(s)*np.cos(gamma)
    c = np.cos(delta)*np.sin(s)*np.sin(gamma)

    return(a,b,c)

def Get_Time_Info(workbook, number):

   # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)

    # Open the Landsat_Input sheet
    ws = workbook['Landsat_Input']

    # Extract Landsat name, number and amount of thermal bands from excel file
    Name_Landsat_Image = str(ws['B%d' %number].value)
    Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)

    # the path to the MTL file of landsat
    Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' % Name_Landsat_Image)

    # read out the general info out of the MTL file
    year, DOY, hour, minutes, UTM_Zone, Sun_elevation = info_general_metadata(Landsat_meta_fileName) # call definition info_general_metadata

    return(year, DOY, hour, minutes, UTM_Zone, Sun_elevation, Landsat_nr)


def Get_LS_Para_Veg(workbook, number, Example_fileName, year, month, day, path_radiance, Apparent_atmosf_transm, cos_zn, dr, vegt_cover_NDVI_max):

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)

    ws = workbook['Additional_Input']

    # If all additional fields are filled in than do not open the datasets
    if ws['B%d' % number].value is None or ws['C%d' % number].value is None:

        print('-------------------- Open Landsat VIS -----------------------')

        # Open the Landsat_Input sheet
        ws = workbook['Landsat_Input']

        # Extract Landsat name, number and amount of thermal bands from excel file
        Name_Landsat_Image = str(ws['B%d' %number].value)
        Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)
                                                             # temperature: 1 = Band 6 for LS_5 & 7, Band 10 for LS_8 (optional)
        # Define bands used for each Landsat number
        if Landsat_nr == 5 or Landsat_nr == 7:
            Bands = np.array([1, 2, 3, 4, 5, 7, 6])
        elif Landsat_nr == 8:
           Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        # Open MTL landsat and get the correction parameters
        Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' %Name_Landsat_Image)
        Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)
        print('Lmin= ', Lmin)
        print('Lmax= ', Lmax)
        print('k1= ', k1_c)
        print('k2= ', k2_c)

        sensor1 = 'LS%d' %Landsat_nr
        res1 = '30m'
        res2 = '30m'
        res3 = '30m'

        # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
        # for the different Landsat images (L5, L7, or L8)
        ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
        ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
        ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])

        # Open one band - To get the metadata of the landsat images only once (to get the extend)
        src_FileName = os.path.join(input_folder, '%s_B2.TIF' %Name_Landsat_Image)  # before 10!
        ls, band_data, ulx, uly, lrx, lry, x_size_ls, y_size_ls = Get_Extend_Landsat(src_FileName)
        print('Original LANDSAT Image - ')
        print('  Size :', x_size_ls, y_size_ls)
        print('  Upper Left corner x, y: ', ulx, ', ', uly)
        print('  Lower right corner x, y: ', lrx, ', ', lry)

        lsc = RC.reproject_dataset_example(src_FileName, Example_fileName)

        #	Get the extend of the remaining landsat file	after clipping based on the DEM file
        y_size_lsc = lsc.RasterYSize
        x_size_lsc = lsc.RasterXSize
        shape_lsc = [x_size_lsc, y_size_lsc]

        print('--- ')
        print('Cropped LANDSAT Image - ')
        print('  Size :', x_size_lsc, y_size_lsc)
        print('  Upper Left corner x, y: ', ulx, ', ',  uly)
        print('  Lower right corner x, y: ', lrx, ', ', lry)

        # if landsat 5 or 7 is used then first create a mask for removing the no data stripes
        if Landsat_nr == 5 or Landsat_nr == 7:
            src_FileName = os.path.join(input_folder, '%s_B6.TIF' % (Name_Landsat_Image)) #open smallest band
            if not os.path.exists(src_FileName):
                src_FileName = os.path.join(input_folder, '%s_B6_VCID_2.TIF' % (Name_Landsat_Image))
            src_FileName_2 = os.path.join(input_folder, '%s_B1.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_3 = os.path.join(input_folder, '%s_B3.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_4 = os.path.join(input_folder, '%s_B4.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_5 = os.path.join(input_folder, '%s_B7.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_6 = os.path.join(input_folder, '%s_B2.TIF' % (Name_Landsat_Image)) #open smallest band
            src_FileName_7 = os.path.join(input_folder, '%s_B5.TIF' % (Name_Landsat_Image)) #open smallest band
            ls_data=Open_landsat(src_FileName,Example_fileName)
            ls_data_2=Open_landsat(src_FileName_2,Example_fileName)
            ls_data_3=Open_landsat(src_FileName_3,Example_fileName)
            ls_data_4=Open_landsat(src_FileName_4,Example_fileName)
            ls_data_5=Open_landsat(src_FileName_5,Example_fileName)
            ls_data_6=Open_landsat(src_FileName_6,Example_fileName)
            ls_data_7=Open_landsat(src_FileName_7,Example_fileName)

            # create and save the landsat mask for all images based on band 11 (smallest map)
            QC_Map=np.zeros((shape_lsc[1], shape_lsc[0]))
            QC_Map=np.where(np.logical_or.reduce((ls_data==0,ls_data_2==0,ls_data_3==0,ls_data_4==0,ls_data_5==0,ls_data_6==0,ls_data_7==0)),1,0)

        # If landsat 8 then use landsat band 10 and 11
        elif Landsat_nr == 8:
             src_FileName_11 = os.path.join(input_folder, '%s_B11.TIF' % (Name_Landsat_Image)) #open smallest band
             ls_data_11=Open_landsat(src_FileName_11,Example_fileName)

             src_FileName_10 = os.path.join(input_folder, '%s_B10.TIF' % (Name_Landsat_Image)) #open smallest band
             ls_data_10=Open_landsat(src_FileName_10, Example_fileName)

             # create and save the landsat mask for all images based on band 10 and 11
             QC_Map=np.zeros((shape_lsc[1], shape_lsc[0]))
             QC_Map=np.where(np.logical_or(ls_data_11==0, ls_data_10==0),1,0)

        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        # Open data of the landsat mask
        ls_data=Open_landsat(src_FileName, Example_fileName)

        # Create 3D array to store Spectral radiance and Reflectivity for each band
        Reflect, Spec_Rad = Landsat_Reflect(Bands, input_folder, Name_Landsat_Image, shape_lsc, QC_Map, Lmax, Lmin, ESUN_L5, ESUN_L7, ESUN_L8, cos_zn, dr, Landsat_nr, Example_fileName)

        # save spectral data
        for i in range(0,6):
            spec_ref_fileName = os.path.join(output_folder, 'Output_radiation_balance','%s_spectral_reflectance_B%s_%s_%s%02d%02d.tif' %(sensor1, Bands[i], res3, year, month, day))
            sebal.save_GeoTiff_proy(lsc, Reflect[:, :, i], spec_ref_fileName, shape_lsc, nband=1)

    else:
        # Get General information example file
        lsc = gdal.Open(Example_fileName)
        nrow = lsc.RasterYSize
        ncol = lsc.RasterXSize
        shape_lsc = [ncol, nrow]

    ######################### Calculate Vegetation Parameters Based on VIS data #####################################

    # Open the Additional input excel sheet
    ws = workbook['Additional_Input']

    # Check NDVI and Calculate NDVI
    try:
        if (ws['B%d' % number].value) is not None:

            # Output folder NDVI
            ndvi_fileName_user = os.path.join(output_folder, 'Output_vegetation', 'User_NDVI_%s_%s%02d%02d.tif' %(res3, year, month, day))
            NDVI=sebal.Reshape_Reproject_Input_data(r'%s' %str(ws['B%d' % number].value),ndvi_fileName_user,Example_fileName)

            water_mask_temp = np.zeros((shape_lsc[1], shape_lsc[0]))
            water_mask_temp[NDVI < 0.0] = 1.0
            sebal.save_GeoTiff_proy(lsc, NDVI, ndvi_fileName_user, shape_lsc, nband=1)

        else:
            # use the Landsat reflectance to calculate the surface albede, NDVI
            NDVI = Calc_NDVI(Reflect)

            # Calculate temporal water mask
            water_mask_temp=sebal.Water_Mask(shape_lsc,Reflect)

    except:
        assert "Please check the NDVI input path"

    # Check Water Mask and replace if it is filled in the additianal data sheet
    try:
        if (ws['E%d' % number].value) is not None:

            # Overwrite the Water mask and change the output name
            water_mask_temp_fileName = os.path.join(output_folder, 'Output_soil_moisture', 'User_Water_mask_temporary_%s_%s%02d%02d.tif' %(res2, year, month, day))
            water_mask_temp = sebal.Reshape_Reproject_Input_data(r'%s' %str(ws['E%d' % number].value), water_mask_temp_fileName, Example_fileName)
            sebal.save_GeoTiff_proy(lsc, water_mask_temp, water_mask_temp_fileName, shape_lsc, nband=1)

    except:
        assert "Please check the Water Mask input path"

    # Check Surface albedo
    try:
        if (ws['C%d' % number].value) is not None:

            # Output folder surface albedo
            surface_albedo_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_albedo_%s_%s%02d%02d.tif' %(res2, year, month, day))
            Surf_albedo=sebal.Reshape_Reproject_Input_data(r'%s' %str(ws['C%d' % number].value),surface_albedo_fileName,Example_fileName)
            sebal.save_GeoTiff_proy(lsc, Surf_albedo, surface_albedo_fileName, shape_lsc, nband=1)

        else:

            # use the Landsat reflectance to calculate the surface albede, NDVI
            Surf_albedo = sebal.Calc_albedo(Reflect, path_radiance, Apparent_atmosf_transm, Landsat_nr)

    except:
          assert "Please check the Albedo input path"

    reflectance_B12 = Reflect[:,:,5]

    # calculate vegetation properties
    FPAR,tir_emis,Nitrogen,vegt_cover,LAI,b10_emissivity=Calc_vegt_para(NDVI, water_mask_temp, vegt_cover_NDVI_max)

    print('Average NDVI = %s' %np.nanmean(NDVI))
    print('Average Surface Albedo = %s' %np.nanmean(Surf_albedo))
    print('Average LAI = %s' %np.nanmean(LAI))
    print('Average Vegetation Cover = %s' %np.nanmean(vegt_cover))
    print('Average FPAR = %s' %np.nanmean(FPAR))

    if not "QC_Map" in locals():
       QC_Map=np.zeros((shape_lsc[1], shape_lsc[0])) 

    return(Surf_albedo, NDVI, LAI, vegt_cover, FPAR, Nitrogen, tir_emis, b10_emissivity, water_mask_temp, QC_Map, reflectance_B12)

def Get_LS_Para_Thermal(workbook, number, Example_fileName, year, month, day, water_mask_temp, b10_emissivity, Temp_inst, Rp, tau_sky, surf_temp_offset, Thermal_Sharpening_not_needed, DEM_fileName, UTM_Zone, eact_inst, QC_Map):

    # import IrriEngine.pySEBAL.pySEBAL_code as SEBAL

    # Open the General input sheet
    ws = workbook['General_Input']

    # Extract the input and output folder, and Image type from the excel file
    input_folder = r"%s" %str(ws['B%d' %number].value)
    output_folder = r"%s" %str(ws['C%d' %number].value)
    Image_Type = 1

    # Open General information example file
    lsc = gdal.Open(Example_fileName)
    nrow = lsc.RasterYSize
    ncol = lsc.RasterXSize
    shape_lsc = [ncol, nrow]

    # Open the Landsat_Input sheet
    ws = workbook['Landsat_Input']

    # Extract Landsat name, number and amount of thermal bands from excel file
    Name_Landsat_Image = str(ws['B%d' %number].value)
    Bands_thermal = int(ws['D%d' %number].value)
    Landsat_nr = int(ws['C%d' %number].value)            # Type of Landsat (LS) image used (LS5, LS7, or LS8)
                                                         # temperature: 1 = Band 6 for LS_5 & 7, Band 10 for LS_8 (optional)
    # Define bands used for each Landsat number
    if Landsat_nr == 5 or Landsat_nr == 7:
        Bands = np.array([1, 2, 3, 4, 5, 7, 6])
    elif Landsat_nr == 8:
        Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
    else:
        print('Landsat image not supported, use Landsat 7 or 8')

    # Open MTL landsat and get the correction parameters
    Landsat_meta_fileName = os.path.join(input_folder, '%s_MTL.txt' %Name_Landsat_Image)
    Lmin, Lmax, k1_c, k2_c = info_band_metadata(Landsat_meta_fileName, Bands)

    sensor1 = 'LS%d' %Landsat_nr
    sensor2 = 'LS%d' %Landsat_nr
    res1 = '30m'
    res2 = '30m'
    res3 = '30m'

    # Open the Landsat_Input sheet
    ws = workbook['Additional_Input']

    # If all additional fields are filled in than do not open the datasets
    if ws['D%d' % number].value is None:

        # Define bands used for each Landsat number
        if Landsat_nr == 5 or Landsat_nr == 7:
            Bands = np.array([1, 2, 3, 4, 5, 7, 6])
        elif Landsat_nr == 8:
           Bands = np.array([2, 3, 4, 5, 6, 7, 10, 11])
        else:
            print('Landsat image not supported, use Landsat 7 or 8')

        print('...................... Open Landsat Thermal ........................')

        # Check if a surface temperature dataset is defined. If so use this one instead of the Landsat, otherwise Landsat
        therm_data = Landsat_therm_data(Bands, input_folder, Name_Landsat_Image, shape_lsc, QC_Map, Example_fileName)

        # Create Cloud mask if BQA map is available (newer version Landsat images)
        BQA_LS_Available = 0
        if os.path.exists(os.path.join(input_folder, '%s_BQA.TIF' %Name_Landsat_Image)):
            src_FileName_BQA = os.path.join(input_folder, '%s_BQA.TIF' %Name_Landsat_Image)
            ls_data_BQA = Open_landsat(src_FileName_BQA, Example_fileName)
            if Landsat_nr == 8:
                Cloud_Treshold = 3700  #2720
            if Landsat_nr == 5 or Landsat_nr == 7:
                Cloud_Treshold = 700
            QC_mask_Cloud = np.copy(ls_data_BQA)
            QC_mask_Cloud[ls_data_BQA<Cloud_Treshold] = 0
            QC_mask_Cloud[ls_data_BQA>=Cloud_Treshold] = 1
            BQA_LS_Available = 1

        # Calculate surface temperature and create a cloud mask
        Surface_temp, cloud_mask_temp = Calc_surface_water_temp(Temp_inst, Landsat_nr, Lmax, Lmin, therm_data, b10_emissivity, k1_c, k2_c, eact_inst, shape_lsc, water_mask_temp, Bands_thermal, Rp, tau_sky, surf_temp_offset, Image_Type)

        # Replace clouds mask calculated by SEBAL by the official BQA file if this exists
        if BQA_LS_Available == 1:
            cloud_mask_temp = QC_mask_Cloud

        Surface_temp[cloud_mask_temp == 1] = np.nan
        print('Mean Surface Temperature = %s Kelvin' %np.nanmean(Surface_temp))

    else:
        try:
            # Output folder surface temperature
            surf_temp_fileName = os.path.join(output_folder, 'Output_vegetation','User_surface_temp_%s_%s%02d%02d.tif' %(res2, year, month, day))
            Surface_temp = sebal.Reshape_Reproject_Input_data(r'%s' %str(ws['D%d' % number].value),surf_temp_fileName,Example_fileName)
            cloud_mask_temp = np.zeros([int(np.shape(Surface_temp)[0]),int(np.shape(Surface_temp)[1])])
            Thermal_Sharpening_not_needed = 0

        except:
            assert "Please check the surface temperature input path"

    return(Surface_temp, cloud_mask_temp, Thermal_Sharpening_not_needed)

#------missing functions
# def Reshape_Reproject_Input_data():
#     return

# def Calc_surface_water_temp():
#     return 

# def save_GeoTiff_proy():
#     return

# def Water_Mask():
#     return

# def Calc_albedo():
#     return

#------------------------------------------------------------------------------
def info_band_metadata(filename, Bands):
    """
    This function retrieves Landsat band information (minimum and maximum
    radiance) from the metadata file.

    """
    Lmin = np.zeros(len(Bands))  # Minimum band radiance, for each band
    Lmax = np.zeros(len(Bands))  # Maximum band radiance, for each band
    k1_const = np.zeros(len(Bands)-6)  # TIRS_Thermal constant k1 ######
    k2_const = np.zeros(len(Bands)-6)  # TIRS_Thermal constant k2 ######
    for band in Bands:
        Landsat_meta = open(filename, "r")  # Open metadata file
        for line in Landsat_meta:
            if re.match("(.*)RADIANCE_MINIMUM_BAND_%1d(.*)" % band, line):
                words = line.split()
                value = float(words[2])
                Lmin[np.where(Bands == band)[0][0]] = value
            if re.match("(.*)RADIANCE_MAXIMUM_BAND_%1d(.*)" % band, line):
                words = line.split()
                value = float(words[2])
                Lmax[np.where(Bands == band)[0][0]] = value
            if re.match("(.*)K1_CONSTANT_BAND_%1d(.*)" % band, line):  # #####
                words = line.split()
                value = float(words[2])
                k1_const[np.where(Bands == band)[0][0]-6] = value
            if re.match("(.*)K2_CONSTANT_BAND_%1d(.*)" % band, line):  # #####
                words = line.split()
                value = float(words[2])
                k2_const[np.where(Bands == band)[0][0]-6] = value
    return Lmin, Lmax, k1_const, k2_const

#------------------------------------------------------------------------------
def Get_Extend_Landsat(src_FileName):
    """
    This function gets the extend of the landsat image
    """
    ls = gdal.Open(src_FileName)       # Open Landsat image
    geo_t_ls = ls.GetGeoTransform()    # Get the Geotransform vector
    x_size_ls = ls.RasterXSize         # Raster xsize - Columns
    y_size_ls = ls.RasterYSize         # Raster ysize - Rows
    (ulx, uly) = geo_t_ls[0], geo_t_ls[3]
    (lrx, lry) = (geo_t_ls[0] + geo_t_ls[1] * x_size_ls,
                  geo_t_ls[3] + geo_t_ls[5] * y_size_ls)
    band_data = ls.GetRasterBand(1)

    return(ls,band_data,ulx,uly,lrx,lry,x_size_ls,y_size_ls)

#------------------------------------------------------------------------------
def Landsat_Reflect(Bands,input_folder,Name_Landsat_Image,shape_lsc,QC_Map,Lmax,Lmin,ESUN_L5,ESUN_L7,ESUN_L8,cos_zn,dr,Landsat_nr, proyDEM_fileName):
    """
    This function calculates and returns the reflectance and spectral radiation from the landsat image.
    """

    Spec_Rad = np.zeros((shape_lsc[1], shape_lsc[0], 7))
    Reflect = np.zeros((shape_lsc[1], shape_lsc[0], 7))
    for band in Bands[:-(len(Bands)-6)]:
        # Open original Landsat image for the band number
        src_FileName = os.path.join(input_folder, '%s_B%1d.TIF'
                                    % (Name_Landsat_Image, band))

        ls_data=Open_landsat(src_FileName, proyDEM_fileName)
        ls_data = ls_data * QC_Map
        # stats = band_data.GetStatistics(0, 1)

        index = np.where(Bands[:-(len(Bands)-6)] == band)[0][0]
        if Landsat_nr == 8:
            # Spectral radiance for each band:
            L_lambda = Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L8, index, cos_zn, dr)
        elif Landsat_nr == 7:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda = Landsat_rho_lambda(L_lambda, ESUN_L7, index, cos_zn, dr)
        elif Landsat_nr == 5:
            # Spectral radiance for each band:
            L_lambda=Landsat_L_lambda(Lmin, Lmax, ls_data, index, Landsat_nr)
            # Reflectivity for each band:
            rho_lambda =Landsat_rho_lambda(L_lambda, ESUN_L5, index, cos_zn, dr)
        else:
            print('Landsat image not supported, use Landsat 5, 7 or 8')

        Spec_Rad[:, :, index] = L_lambda
        Reflect[:, :, index] = rho_lambda
    Reflect = Reflect.clip(0.0, 1.0)
    return(Reflect,Spec_Rad)


#------------------------------------------------------------------------------
def Landsat_L_lambda(Lmin,Lmax,ls_data,index,Landsat_nr):
    """
    Calculates the lambda from landsat
    """
    if Landsat_nr==8:
        L_lambda = ((Lmax[index] - Lmin[index]) / (65535 - 1) * ls_data + Lmin[index])
    elif Landsat_nr == 5 or Landsat_nr ==7:
        L_lambda = (Lmax[index] - Lmin[index]) / 255 * ls_data + Lmin[index]
    return(L_lambda)


#------------------------------------------------------------------------------
def Landsat_rho_lambda(L_lambda,ESUN,index,cos_zn,dr):
    """
    Calculates the rho from landsat
    """
    rho_lambda = np.pi * L_lambda / (ESUN[index] * cos_zn * dr)
    return(rho_lambda)


#------------------------------------------------------------------------------
def Landsat_therm_data(Bands, input_folder, Name_Landsat_Image, shape_lsc, QC_Map, proyDEM_fileName):
    """
    This function calculates and returns the thermal data from the landsat image.
    """

    therm_data = np.zeros((shape_lsc[1], shape_lsc[0], len(Bands)-6))
    for band in Bands[-(len(Bands)-6):]:
        
        # Open original Landsat image for the band number
        src_FileName = os.path.join(input_folder, '%s_B%1d.TIF'
                                    % (Name_Landsat_Image, band))
        if not os.path.exists(src_FileName):
             src_FileName = os.path.join(input_folder, '%s_B%1d_VCID_2.TIF'
                                    % (Name_Landsat_Image, band))

        ls_data = Open_landsat(src_FileName, proyDEM_fileName)
        ls_data = ls_data * QC_Map
        index = np.where(Bands[:] == band)[0][0] - 6
        therm_data[:, :, index] = ls_data

    return(therm_data)

#------------------------------------------------------------------------------
def Open_landsat(src_FileName, proyDEM_fileName):
    """
    This function opens a landsat image and returns the data array of a specific landsat band.
    """

    # crop band to the DEM extent
    ls = RC.reproject_dataset_example(src_FileName, proyDEM_fileName)

    # Open the cropped Landsat image for the band number
    ls_data = ls.GetRasterBand(1).ReadAsArray()
    return(ls_data)

#%%

if __name__ == "__main__":
    LS_input_folder = r"/Volumes/Data/RS_scenes/india_case/C1_L1"
    output_folder = r"/Volumes/Data/RS_scenes/india_case/C1_L1_output"

    # Find dates
    os.chdir(LS_input_folder)
    files_LS = glob.glob("*.tar.gz")
    dates = [datetime.datetime.strptime(k.split("_")[3], "%Y%m%d") for k in files_LS]

    for date in dates:

        date_str = date.strftime("%Y%m%d")
        os.chdir(LS_input_folder)
        found = glob.glob("L*_*_%s_*_*_*.tar.gz" %date_str)
        
        if len(found) > 0:
            
            for found_one in found:

                extract_folder = os.path.join(LS_input_folder, found_one.replace(".tar.gz",""))
                if not os.path.exists(extract_folder):
                    os.makedirs(extract_folder)

                os.chdir(extract_folder)
                
                LS_file = os.path.splitext(os.path.splitext(found_one)[0])[0] +"_MTL.txt"
                Startdate_IrriEngine = datetime.datetime.strptime(LS_file.split("_")[3], "%Y%m%d")
                Startdate_IrriEngine = Startdate_IrriEngine.strftime("%Y-%m-%d")
                    
                file_LS_Cloud = os.path.join(output_folder, "%s" %LS_file.split("_")[2], "%s"  %LS_file.split("_")[3], "LS_CLOUD_%s.tif" %("_".join(LS_file.split("_")[0:4])))
                Tile_name = LS_file.split("_")[2]
                
                if not os.path.exists(file_LS_Cloud):
        
                    print("Process %s" %found_one)

                    LS_B1_filename = LS_file.replace("_MTL.txt", "_B1.tif")

                    if not os.path.exists(LS_B1_filename):
                        # Extract zip file
                        DC.Extract_Data_tar_gz(os.path.join(LS_input_folder,found_one), extract_folder)
                
                    # Get LS_filename
                    print(found_one)
                    Landsat_nr = int(found_one[3])
                    
                    # read out the general info out of the MTL file
                    year, DOY, hour, minutes, UTM_Zone, Sun_elevation = info_general_metadata(LS_file)     
                    UTM_Zone = int("326%02d" %UTM_Zone)#!!! nu alleen voor noord
                    folder_LS_RAW = os.path.join(LS_input_folder, "RAW")  
                    folder_LS_RAW_DEM = os.path.join(folder_LS_RAW, "%s" %(Tile_name))         
                    filename_DEM = os.path.join(folder_LS_RAW_DEM, "HydroSHED", "DEM", "DEM_HydroShed_m_3s.tif")
                        
                    if not os.path.exists(filename_DEM):
                        
                        if not os.path.exists(folder_LS_RAW_DEM):
                            os.makedirs(folder_LS_RAW_DEM)
                        
                        # Get extend in degrees of tile
                        dest = gdal.Open(os.path.join(extract_folder,LS_B1_filename))
                        geo = dest.GetGeoTransform()
                        epsg_from = Proj(init="epsg:%s" %UTM_Zone) 
                        epsg_to = Proj(init="epsg:4326")
                        x1 = geo[0]
                        x2 = geo[0] + dest.RasterXSize*geo[1]
                        y1 = geo[3]
                        y2 = geo[3] + dest.RasterYSize * geo[5]
                        x1_deg, y1_deg = transform(epsg_from, epsg_to, x1, y1)
                        x2_deg, y2_deg = transform(epsg_from, epsg_to, x2, y2)
                        
                        watertools.Collect.DEM.HydroSHED(folder_LS_RAW_DEM, [y2_deg - 0.1, y1_deg + 0.1], [x1_deg - 0.1, x2_deg + 0.1])
                        
                    else:
                            # Get extend in degrees of tile
                        dest = gdal.Open(os.path.join(extract_folder,LS_B1_filename))
                        geo = dest.GetGeoTransform()
                        epsg_from = Proj(init="epsg:%s" %UTM_Zone) 
                        epsg_to = Proj(init="epsg:4326")
                        x1 = geo[0]
                        x2 = geo[0] + dest.RasterXSize*geo[1]
                        y1 = geo[3]
                        y2 = geo[3] + dest.RasterYSize * geo[5]
                        x1_deg, y1_deg = transform(epsg_from, epsg_to, x1, y1)
                        x2_deg, y2_deg = transform(epsg_from, epsg_to, x2, y2)
                        
                    lon_ave = (x1_deg + x2_deg)/2  
                        
                    GMT_offset = round(lon_ave * 24 / 360)
                    
                    # Find Time to take
                    time_wheater_data = np.linspace(0,21,8) + 1.5

                    # Collect required METEO data
                    Time_Step = np.argwhere(time_wheater_data == Find_nearest_time(time_wheater_data, hour + minutes/60))[0][0]
                    hour_take = Time_Step * 3
                    
                    temp_format_inst = os.path.join(folder_LS_RAW, "Weather_Data", "Model", "GLDAS", "three_hourly", "tair_f_inst", "Tair_GLDAS-NOAH_C_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
                    RH_format_inst = os.path.join(folder_LS_RAW, "Weather_Data", "Model", "GLDAS", "three_hourly", "rh_f_inst", "Hum_GLDAS-NOAH_percentage_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
                
                    if not os.path.exists(temp_format_inst.format(yyyy=date.year, mm=date.month, dd=date.day, HH=hour_take)):
                        watertools.Collect.GLDAS.three_hourly(folder_LS_RAW, ['tair_f_inst'], date, date, [y2_deg - 0.3, y1_deg + 0.3], [x1_deg - 0.3, x2_deg + 0.3], Periods = [int(Time_Step+1)])
                                        
                    if not os.path.exists(RH_format_inst.format(yyyy=date.year, mm=date.month, dd=date.day, HH=hour_take)):          
                            watertools.Collect.GLDAS.three_hourly(folder_LS_RAW, ['psurf_f_inst', 'qair_f_inst'], date, date, [y2_deg - 0.3, y1_deg + 0.3], [x1_deg - 0.3, x2_deg + 0.3], Periods = [int(Time_Step+1)])
                            input_format_psurf = os.path.join(folder_LS_RAW, "Weather_Data", "Model", "GLDAS", "three_hourly", "psurf_f_inst", "P_GLDAS-NOAH_kpa_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
                            input_format_qair = os.path.join(folder_LS_RAW, "Weather_Data", "Model", "GLDAS", "three_hourly", "qair_f_inst", "Hum_GLDAS-NOAH_kg-kg_3hour_{yyyy}.{mm:02d}.{dd:02d}_{HH:02d}00.tif")
                            watertools.Products.RH.Use_Meteo.Calc_Humidity(temp_format_inst,
                                                                input_format_psurf,
                                                                input_format_qair,
                                                                RH_format_inst,
                                                                date, date+pd.DateOffset(days=1), freq = "3H")                
                        
                    # Define output folders
                    folder_LS_NDVI = os.path.join(output_folder, "%s" %(Tile_name), "NDVI")
                    folder_LS_Albedo = os.path.join(output_folder, "%s" %(Tile_name), "Albedo")
                    folder_LS_VSDI = os.path.join(output_folder, "%s" %(Tile_name), "VSDI")
                    folder_LS_NDII = os.path.join(output_folder, "%s" %(Tile_name), "NDII")
                    folder_LS_B12 = os.path.join(output_folder, "%s" %(Tile_name), "B12")
                    folder_LS_LST = os.path.join(output_folder, "%s" %(Tile_name), "LST")
                    folder_LS_Cloud = os.path.join(output_folder, "%s" %(Tile_name), "CLOUD")
                    
                    # Save output folders
                    FOLDERS = [extract_folder, folder_LS_RAW, folder_LS_NDVI, folder_LS_Albedo, folder_LS_VSDI, folder_LS_NDII, folder_LS_B12, folder_LS_Cloud, folder_LS_LST]

                    for folder in FOLDERS:
                        if not os.path.exists(folder) and folder != '':
                            print("create input folder %s" %folder)
                            os.makedirs(folder)
                            
                    Process_LS_Image(LS_file, FOLDERS, filename_DEM, UTM_Zone, Startdate_IrriEngine, GMT_offset)

                    shutil.rmtree(extract_folder)

#%%
