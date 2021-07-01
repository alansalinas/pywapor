# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:43:42 2020

@author: timhe
"""

import os
import glob
import gdal
import datetime 
import numpy as np
import pandas as pd
from pyproj import Proj, transform

import watertools
import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC
# import IrriEngine
# import IrriEngine.Constants as Constants
# import IrriEngine.pySEBAL.pySEBAL_input_LANDSAT as SEBAL_LS
# import IrriEngine.pySEBAL.pySEBAL_code as SEBAL
import SEBAL.pySEBAL.pySEBAL_input_LANDSAT as SEBAL_LS
import SEBAL.pySEBAL.pySEBAL_code as SEBAL
class ConstantsObject:
    # between 0.025 and 0.040, average portion of the incoming 
    # solar radiation across all bands that is back-scattered 
    # to the satellite before it reaches the earths surface.
    path_radiance_LS = 0.03
    # between 0.75 and 0.9, atmospheric transmissivity = 0.75 + 2e-5 * z.
    Apparent_atmosf_transm_LS = 0.8 
    # Use this for uncorrected radiance.
    Rp = 0.0 
    tau_sky = 1.0
    # Used for a correction of the cloudmask which is not used.
    surf_temp_offset = -9999
Constants = ConstantsObject()

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
    year, DOY, hour_GTM, minutes_GTM, UTM_Zone_LS, Sun_elevation = SEBAL_LS.info_general_metadata(Landsat_meta_fileName)    

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
                lat_proy_UTM, lon_proy_UTM  = SEBAL.DEM_lat_lon(dest_dem)      
            
                # Calculate Latitude and Longitude files in WGS84
                input_projection = Proj(init="epsg:%s" %UTM_Zone)
                output_projection = Proj(init="epsg:4326")
                lon_proy, lat_proy = transform(input_projection, output_projection, lon_proy_UTM.flatten(), lat_proy_UTM.flatten())
                lon_proy.resize((nrow, ncol))
                lat_proy.resize((nrow, ncol))
     
                # Calculate aspect and slope
                slope, aspect = SEBAL.Calc_Gradient(DEM_resh, 30)
    
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
            Lmin, Lmax, k1_c, k2_c = SEBAL_LS.info_band_metadata(Landsat_meta_fileName, Bands)
    
            # Mean solar exo-atmospheric irradiance for each band (W/m2/microm)
            # for the different Landsat images (L5, L7, or L8)
            ESUN_L5 = np.array([1983, 1796, 1536, 1031, 220, 83.44])
            ESUN_L7 = np.array([1997, 1812, 1533, 1039, 230.8, 84.9])
            ESUN_L8 = np.array([1973.28, 1842.68, 1565.17, 963.69, 245, 82.106])   
            
            # Create MASK for no data values LS
            filename_BQA = os.path.join(FOLDERS[0], '%s_BQA.TIF' % (Name_Landsat_Image))
            LS_BQA = SEBAL_LS.Open_landsat(filename_BQA, dest_dem)
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
            Reflect, Spec_Rad = SEBAL_LS.Landsat_Reflect(Bands, FOLDERS[0], Name_Landsat_Image, shape_dem, QC_mask_Cloud, Lmax, Lmin, ESUN_L5, ESUN_L7, ESUN_L8, cos_zn, dr, Landsat_nr, dest_dem)
            
            # Calculate NDVI
            print("calculate spectral indices")                  
            NDVI = SEBAL.Calc_NDVI(Reflect)
            NDII = (Reflect[:, :, 3] - Reflect[:, :, 4])/(Reflect[:, :, 3] + Reflect[:, :, 4])
            VSDI = 1-((Reflect[:, :, 4] - Reflect[:, :, 0]) + (Reflect[:, :, 2] - Reflect[:, :, 0]))
        
            # Surface albedo:
            Surf_albedo = (0.3 * Reflect[:, :, 0] + 0.277 * Reflect[:, :, 1] +
                           0.233 * Reflect[:, :, 2] + 0.143 * Reflect[:, :, 3] +
                           0.036 * Reflect[:, :, 4] + 0.012 * Reflect[:, :, 5] -
                           Constants.path_radiance_LS) / np.power(Constants.Apparent_atmosf_transm_LS, 2)            
            
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
            FPAR, tir_emis, Nitrogen, vegt_cover, LAI, b10_emissivity = SEBAL.Calc_vegt_para(NDVI, water_mask_temp)
    
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
                therm_data = SEBAL_LS.Landsat_therm_data(Bands, FOLDERS[0], Name_Landsat_Image, shape_dem, QC_mask_Cloud, dest_dem)      
                therm_data[therm_data==0] = np.nan
                Surface_temp, cloud_mask_temp = SEBAL.Calc_surface_water_temp(Temperature_inst, Landsat_nr, Lmax, Lmin, therm_data, b10_emissivity, k1_c, k2_c, eact_inst, shape_dem, water_mask_temp, Bands_thermal, Constants.Rp, Constants.tau_sky, Constants.surf_temp_offset, 1)
          
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
    
def Create_Example_UTM_File(filename_ex, pixel_spacing, proj):

    # Reproject DEM to fit LS scene
    dest, ulx, lry, lrx, uly, epsg_to = SEBAL.reproject_dataset(
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
    w = SEBAL.w_time(GMT_time, lon_proy, DOY)                            # Hour angle (radians)
    a,b,c = SEBAL.Constants(delta,s,gamma,phi)
    cos_zn= SEBAL.AngleSlope(a,b,c,w)
    cos_zn = cos_zn.clip(Min_cos_zn, Max_cos_zn)

    print('Average Cos Zenith Angle: ', '%0.3f (Radians)' % np.nanmean(cos_zn))

    dr = 1 + 0.033 * np.cos(DOY*2*np.pi/365)  # Inverse relative distance Earth-Sun

    return(cos_zn, dr)    

LS_input_folder = r"F:\Project_FAOJORDAN\Input_Data\LS"
output_folder = r"F:\Project_FAOJORDAN\Input_Data\LS_input_WAPOR"

# Find dates
os.chdir(LS_input_folder)
files_LS = glob.glob("*.tar.gz")
dates = [datetime.datetime.strptime(k.split("_")[3], "%Y%m%d") for k in files_LS]

for date in dates:
    
    if date.year == 2015:
        date_str = date.strftime("%Y%m%d")
        os.chdir(LS_input_folder)
        found = glob.glob("L*_*_%s_*_*_*.tar.gz" %date_str)
        
        if len(found) > 0:
            
            for found_one in found:
                
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
                        DC.Extract_Data_tar_gz(os.path.join(LS_input_folder,found_one), LS_input_folder)
             
                    # Get LS_filename
                    print(found_one)
                    Landsat_nr = int(found_one[3])
                    
                    # read out the general info out of the MTL file
                    year, DOY, hour, minutes, UTM_Zone, Sun_elevation = SEBAL_LS.info_general_metadata(LS_file)     
                    UTM_Zone = int("326%02d" %UTM_Zone)#!!! nu alleen voor noord
                    folder_LS_RAW = os.path.join(LS_input_folder, "RAW")  
                    folder_LS_RAW_DEM = os.path.join(folder_LS_RAW, "%s" %(Tile_name))         
                    filename_DEM = os.path.join(folder_LS_RAW_DEM, "HydroSHED", "DEM", "DEM_HydroShed_m_3s.tif")
                       
                    if not os.path.exists(filename_DEM):
                        
                        if not os.path.exists(folder_LS_RAW_DEM):
                            os.makedirs(folder_LS_RAW_DEM)
                        
                        # Get extend in degrees of tile
                        dest = gdal.Open(os.path.join(LS_input_folder,LS_B1_filename))
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
                        dest = gdal.Open(os.path.join(LS_input_folder,LS_B1_filename))
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
                    time= hour + minutes/60
                    time_wheater_data = np.linspace(0,21,8) + 1.5
    
                    # Collect required METEO data
                    Time_Step = np.argwhere(time_wheater_data == Find_nearest_time(time_wheater_data, time))[0][0]
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
                    FOLDERS = [LS_input_folder, folder_LS_RAW, folder_LS_NDVI, folder_LS_Albedo, folder_LS_VSDI, folder_LS_NDII, folder_LS_B12, folder_LS_Cloud, folder_LS_LST]
    
                    for folder in FOLDERS:
                        if not os.path.exists(folder) and folder != '':
                            print("create input folder %s" %folder)
                            os.makedirs(folder)
                            
                    Process_LS_Image(LS_file, FOLDERS, filename_DEM, UTM_Zone, Startdate_IrriEngine, GMT_offset)