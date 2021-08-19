# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:31:41 2020

@author: timhe
"""


import os
from osgeo import gdal
import datetime
import numpy as np

import watertools.General.data_conversions as DC
import watertools.General.raster_conversions as RC


#LS_folder = r"F:\Project_FAOJORDAN\Input_Data\LS_output_WAPOR"
#composite_folder = r"D:\Project_ISRPAL\SWEO_Input\Variables\ETLook\ETLook_output_composite"
#output_folder = r"F:\Project_FAOJORDAN\Input_Data\LS_output_WAPOR_composite"

#LS_LVL1_maps = os.listdir(composite_folder)
#DATES_LS_LVL1 = np.array([datetime.datetime.strptime(k, "%Y%m%d").toordinal() for k in LS_LVL1_maps if len(k) == 8])

#for DATE_LS_LVL1 in DATES_LS_LVL1:
#    Date = datetime.datetime.fromordinal(DATE_LS_LVL1)
#    if Date.year==2015:
#       main_et(Date, LS_folder, composite_folder, output_folder)
       #main_e(Date, LS_folder, composite_folder, output_folder)
       #main_bio(Date, LS_folder, composite_folder, output_folder)
       #main_sm(Date, LS_folder, composite_folder, output_folder)       
       
def main_et(Date, LS_folder, composite_folder, output_folder):

    date_comp = Date.strftime("%Y%m%d")
    filename_end_out = os.path.join(output_folder, "et_24_mm_%s.tif" %(date_comp))   
 
    if not os.path.exists(filename_end_out):
        
        # Find all LS dates
        LS_LVL3_maps = os.listdir(LS_folder)
        DATES_LS_LVL3 = np.array([datetime.datetime.strptime(k, "%Y%m%d").toordinal() for k in LS_LVL3_maps if len(k) == 8])
        
        # Find amount of days between
        Diff_days = np.abs(DATES_LS_LVL3 - Date.toordinal())
        Diff_days_sorted = np.sort(Diff_days)
        
        # Create LST_rel map
        output_folder_LSTrel = os.path.join(LS_folder, "LSTrel")
        if not os.path.exists(output_folder_LSTrel):
            os.makedirs(output_folder_LSTrel)
            
        # 
        Continue = 1
        
        for Diff_day_sorted in [k for k in Diff_days_sorted if k <60]:
            
            if Continue == 1:
                
                take_days = np.argwhere(Diff_days == Diff_day_sorted)
                
                for take_day in take_days:
                
                    DATE_LS_LVL3 = int(DATES_LS_LVL3[take_day])
                    
                    date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")
                
                    print("Add %s in relative map" %date_LS)
                    
                    Output_filename_LSTrel = os.path.join(output_folder_LSTrel, "lstrel_%s.tif" %date_LS)
                    
                    if not os.path.exists(Output_filename_LSTrel):
                    
                        ET_map_LS = os.path.join(LS_folder, date_LS, "et_24_mm_%s.tif" %(date_LS))
                        
                        dest_LS_30 = gdal.Open(ET_map_LS)
                        geo_30 = dest_LS_30.GetGeoTransform()
                        proj = dest_LS_30.GetProjection()
                        size_x = dest_LS_30.RasterXSize
                        size_y = dest_LS_30.RasterYSize
                        Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
                        
                        geo_1000 = np.array(geo_30)
                        geo_1000[1] = geo_1000[1]*108
                        geo_1000[5] = geo_1000[5]*108
                        geo_1000 = tuple(geo_1000)
                        
                        # Create MEM up file
                        dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
                        dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)
                        
                        # upscale LS
                        # dest_LS_up = RC.reproject_dataset_example(dest_LS_30, dest_up, 4)
            
                        # Get relative array
                        relMap_LS_one = Create_Relative_Map(Array_LS, 108, Include_NAN = 0)  
                        
                        DC.Save_as_tiff(Output_filename_LSTrel, relMap_LS_one, geo_30, proj)
                    
                    else:
                        dest_lstrel = gdal.Open(Output_filename_LSTrel)
                        relMap_LS_one = dest_lstrel.GetRasterBand(1).ReadAsArray()
                    
                    if "relMap_LS" not in locals():
                        relMap_LS = np.ones(relMap_LS_one.shape) * np.nan
                        
                    # Fit in endmap
                    relMap_LS = np.where(np.isnan(relMap_LS), relMap_LS_one, relMap_LS)
                
                    Amount_Of_Pixels2Go = np.count_nonzero(np.isnan(relMap_LS))
                
                if Amount_Of_Pixels2Go == 0:
                    Continue = 0
                else:
                    print("Still %s pixels to go" %Amount_Of_Pixels2Go)

        if "dest_up" not in locals():
            
            DATE_LS_LVL3 = int(DATES_LS_LVL3[0])   
            date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")           
            ET_map_LS = os.path.join(LS_folder, date_LS, "et_24_mm_%s.tif" %(date_LS))
                        
            dest_LS_30 = gdal.Open(ET_map_LS)
            geo_30 = dest_LS_30.GetGeoTransform()
            proj = dest_LS_30.GetProjection()
            size_x = dest_LS_30.RasterXSize
            size_y = dest_LS_30.RasterYSize
            Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
            
            geo_1000 = np.array(geo_30)
            geo_1000[1] = geo_1000[1]*108
            geo_1000[5] = geo_1000[5]*108
            geo_1000 = tuple(geo_1000)
            
            # Create MEM up file
            dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
            dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)   
         
        ET_map_Comp = os.path.join(composite_folder, date_comp, "et_24_mm_%s.tif" %(date_comp)) 
        dest_comp_up = RC.reproject_dataset_example(ET_map_Comp, dest_up, 4) 
        
        # Downscale back to 30m
        dest_comp_down = RC.reproject_dataset_example(dest_comp_up, dest_LS_30, 2) 
        Array_comp = dest_comp_down.GetRasterBand(1).ReadAsArray()
                                  
        Temp_array = relMap_LS.clip(0.0001, 6) * Array_comp
        
        # Scale Temp_array
        dest_mem_temp = DC.Save_as_MEM(Temp_array, geo_30, 4326)
        
        # Open input composite and get geotransform
        dest_comp = gdal.Open(ET_map_Comp)
        geo_comp = dest_comp.GetGeoTransform()
        proj_comp = dest_comp.GetProjection()
        Array_comp_MODIS = dest_comp.GetRasterBand(1).ReadAsArray()
        Array_comp_MODIS_BOX = Create_Relative_Map(Array_comp_MODIS, 3, Include_NAN = 1, type_return = "mean")
        
        # Upscale LS comp
        dest_LScomp_up = RC.reproject_dataset_example(dest_mem_temp, ET_map_Comp, 4)
        Array_LScomp_up = dest_LScomp_up.GetRasterBand(1).ReadAsArray()
        Array_LScomp_up_BOX = Create_Relative_Map(Array_LScomp_up, 3, Include_NAN = 1, type_return = "mean")
        
        # Calculate Fraction Correction
        Corr_Frac = Array_comp_MODIS_BOX/Array_LScomp_up_BOX
        Corr_Frac = Corr_Frac.clip(0.0001, 6)
        dest_corr_up = DC.Save_as_MEM(Corr_Frac, geo_comp, proj_comp)
        
        # downscale corr factor
        dest_corr_down = RC.reproject_dataset_example(dest_corr_up, dest_mem_temp, 2)
        Corr_down = dest_corr_down.GetRasterBand(1).ReadAsArray()
        
        # apply corr
        End_array = Corr_down * Temp_array
        End_array[End_array>11] = np.nan
        
        DC.Save_as_tiff(filename_end_out, End_array, geo_30, 4326)
    
    return()
                
  
def main_e(Date, LS_folder, composite_folder, output_folder):

    date_comp = Date.strftime("%Y%m%d")
    filename_end_out = os.path.join(output_folder, "e_24_mm_%s.tif" %(date_comp))   
 
    if not os.path.exists(filename_end_out):
        
        # Find all LS dates
        LS_LVL3_maps = os.listdir(LS_folder)
        DATES_LS_LVL3 = np.array([datetime.datetime.strptime(k, "%Y%m%d").toordinal() for k in LS_LVL3_maps if len(k) == 8])
        
        # Find amount of days between
        Diff_days = np.abs(DATES_LS_LVL3 - Date.toordinal())
        Diff_days_sorted = np.sort(Diff_days)
        
        # Create E_rel map
        output_folder_Erel = os.path.join(LS_folder, "Erel")
        if not os.path.exists(output_folder_Erel):
            os.makedirs(output_folder_Erel)
            
        # 
        Continue = 1
        
        for Diff_day_sorted in [k for k in Diff_days_sorted if k <60]:
            
            if Continue == 1:
                
                take_days = np.argwhere(Diff_days == Diff_day_sorted)
                
                for take_day in take_days:
                
                    DATE_LS_LVL3 = int(DATES_LS_LVL3[take_day])
                    
                    date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")
                
                    print("Add %s in relative map" %date_LS)
                    
                    Output_filename_Erel = os.path.join(output_folder_Erel, "erel_%s.tif" %date_LS)
                    
                    if not os.path.exists(Output_filename_Erel):
                    
                        E_map_LS = os.path.join(LS_folder, date_LS, "e_24_mm_%s.tif" %(date_LS))
                        
                        dest_LS_30 = gdal.Open(E_map_LS)
                        geo_30 = dest_LS_30.GetGeoTransform()
                        proj = dest_LS_30.GetProjection()
                        size_x = dest_LS_30.RasterXSize
                        size_y = dest_LS_30.RasterYSize
                        Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
                        
                        geo_1000 = np.array(geo_30)
                        geo_1000[1] = geo_1000[1]*108
                        geo_1000[5] = geo_1000[5]*108
                        geo_1000 = tuple(geo_1000)
                        
                        # Create MEM up file
                        dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
                        dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)
                        
                        # upscale LS
                        dest_LS_up = RC.reproject_dataset_example(dest_LS_30, dest_up, 4)
            
                        # Get relative array
                        relMap_LS_one = Create_Relative_Map(Array_LS, 108, Include_NAN = 0)  
                        
                        DC.Save_as_tiff(Output_filename_Erel, relMap_LS_one, geo_30, proj)
                    
                    else:
                        dest_erel = gdal.Open(Output_filename_Erel)
                        relMap_LS_one = dest_erel.GetRasterBand(1).ReadAsArray()
                    
                    if "relMap_LS" not in locals():
                        relMap_LS = np.ones(relMap_LS_one.shape) * np.nan
                        
                    # Fit in endmap
                    relMap_LS = np.where(np.isnan(relMap_LS), relMap_LS_one, relMap_LS)
                
                    Amount_Of_Pixels2Go = np.count_nonzero(np.isnan(relMap_LS))
                
                if Amount_Of_Pixels2Go == 0:
                    Continue = 0
                else:
                    print("Still %s pixels to go" %Amount_Of_Pixels2Go)
         
        if "dest_up" not in locals():
            
            DATE_LS_LVL3 = int(DATES_LS_LVL3[0])   
            date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")           
            E_map_LS = os.path.join(LS_folder, date_LS, "e_24_mm_%s.tif" %(date_LS))
                        
            dest_LS_30 = gdal.Open(E_map_LS)
            geo_30 = dest_LS_30.GetGeoTransform()
            proj = dest_LS_30.GetProjection()
            size_x = dest_LS_30.RasterXSize
            size_y = dest_LS_30.RasterYSize
            Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
            
            geo_1000 = np.array(geo_30)
            geo_1000[1] = geo_1000[1]*108
            geo_1000[5] = geo_1000[5]*108
            geo_1000 = tuple(geo_1000)
            
            # Create MEM up file
            dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
            dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)   
         
        E_map_Comp = os.path.join(composite_folder, date_comp, "e_24_mm_%s.tif" %(date_comp)) 
        dest_comp_up = RC.reproject_dataset_example(E_map_Comp, dest_up, 4) 
        
        # Downscale back to 30m
        dest_comp_down = RC.reproject_dataset_example(dest_comp_up, dest_LS_30, 2) 
        Array_comp = dest_comp_down.GetRasterBand(1).ReadAsArray()
                                  
        Temp_array = relMap_LS.clip(0.0001, 6) * Array_comp
        
        # Scale Temp_array
        dest_mem_temp = DC.Save_as_MEM(Temp_array, geo_30, 4326)
        
        # Open input composite and get geotransform
        dest_comp = gdal.Open(E_map_Comp)
        geo_comp = dest_comp.GetGeoTransform()
        proj_comp = dest_comp.GetProjection()
        Array_comp_MODIS = dest_comp.GetRasterBand(1).ReadAsArray()
        Array_comp_MODIS_BOX = Create_Relative_Map(Array_comp_MODIS, 3, Include_NAN = 1, type_return = "mean")
        
        # Upscale LS comp
        dest_LScomp_up = RC.reproject_dataset_example(dest_mem_temp, E_map_Comp, 4)
        Array_LScomp_up = dest_LScomp_up.GetRasterBand(1).ReadAsArray()
        Array_LScomp_up_BOX = Create_Relative_Map(Array_LScomp_up, 3, Include_NAN = 1, type_return = "mean")
        
        # Calculate Fraction Correction
        Corr_Frac = Array_comp_MODIS_BOX/Array_LScomp_up_BOX
        Corr_Frac = Corr_Frac.clip(0.0001, 6)
        dest_corr_up = DC.Save_as_MEM(Corr_Frac, geo_comp, proj_comp)
        
        # downscale corr factor
        dest_corr_down = RC.reproject_dataset_example(dest_corr_up, dest_mem_temp, 2)
        Corr_down = dest_corr_down.GetRasterBand(1).ReadAsArray()
        
        # apply corr
        End_array = Corr_down * Temp_array
        End_array[End_array>11] = np.nan
        
        DC.Save_as_tiff(filename_end_out, End_array, geo_30, 4326)
    
    return()
                

            
def main_bio(Date, LS_folder, composite_folder, output_folder):

    date_comp = Date.strftime("%Y%m%d")
    filename_end_out = os.path.join(output_folder, "biomass_prod_kg-ha_%s.tif" %(date_comp))   
 
    if not os.path.exists(filename_end_out):
        
        # Find all LS dates
        LS_LVL3_maps = os.listdir(LS_folder)
        DATES_LS_LVL3 = np.array([datetime.datetime.strptime(k, "%Y%m%d").toordinal() for k in LS_LVL3_maps if len(k) == 8])
        
        # Find amount of days between
        Diff_days = np.abs(DATES_LS_LVL3 - Date.toordinal())
        Diff_days_sorted = np.sort(Diff_days)
        
        # Create BIO_rel map
        output_folder_BIOrel = os.path.join(LS_folder, "BIOrel")
        if not os.path.exists(output_folder_BIOrel):
            os.makedirs(output_folder_BIOrel)
            
        # 
        Continue = 1
        
        for Diff_day_sorted in [k for k in Diff_days_sorted if k <60]:
            
            if Continue == 1:
                
                take_days = np.argwhere(Diff_days == Diff_day_sorted)
                
                for take_day in take_days:
                
                    DATE_LS_LVL3 = int(DATES_LS_LVL3[take_day])
                    
                    date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")
                
                    print("Add %s in relative map" %date_LS)
                    
                    Output_filename_BIOrel = os.path.join(output_folder_BIOrel, "biorel_%s.tif" %date_LS)
                    
                    if not os.path.exists(Output_filename_BIOrel):
                    
                        BIO_map_LS = os.path.join(LS_folder, date_LS, "biomass_prod_kg-ha_%s.tif" %(date_LS))
                        
                        dest_LS_30 = gdal.Open(BIO_map_LS)
                        geo_30 = dest_LS_30.GetGeoTransform()
                        proj = dest_LS_30.GetProjection()
                        size_x = dest_LS_30.RasterXSize
                        size_y = dest_LS_30.RasterYSize
                        Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
                        
                        geo_1000 = np.array(geo_30)
                        geo_1000[1] = geo_1000[1]*108
                        geo_1000[5] = geo_1000[5]*108
                        geo_1000 = tuple(geo_1000)
                        
                        # Create MEM up file
                        dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
                        dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)
                        
                        # upscale LS
                        dest_LS_up = RC.reproject_dataset_example(dest_LS_30, dest_up, 4)
            
                        # Get relative array
                        relMap_LS_one = Create_Relative_Map(Array_LS, 108, Include_NAN = 0)  
                        
                        DC.Save_as_tiff(Output_filename_BIOrel, relMap_LS_one, geo_30, proj)
                    
                    else:
                        dest_biorel = gdal.Open(Output_filename_BIOrel)
                        relMap_LS_one = dest_biorel.GetRasterBand(1).ReadAsArray()
                    
                    if "relMap_LS" not in locals():
                        relMap_LS = np.ones(relMap_LS_one.shape) * np.nan
                        
                    # Fit in endmap
                    relMap_LS = np.where(np.isnan(relMap_LS), relMap_LS_one, relMap_LS)
                
                    Amount_Of_Pixels2Go = np.count_nonzero(np.isnan(relMap_LS))
                
                if Amount_Of_Pixels2Go == 0:
                    Continue = 0
                else:
                    print("Still %s pixels to go" %Amount_Of_Pixels2Go)
         
        if "dest_up" not in locals():
            
            DATE_LS_LVL3 = int(DATES_LS_LVL3[0])   
            date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")           
            BIO_map_LS = os.path.join(LS_folder, date_LS, "biomass_prod_kg-ha_%s.tif" %(date_LS))
                        
            dest_LS_30 = gdal.Open(BIO_map_LS)
            geo_30 = dest_LS_30.GetGeoTransform()
            proj = dest_LS_30.GetProjection()
            size_x = dest_LS_30.RasterXSize
            size_y = dest_LS_30.RasterYSize
            Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
            
            geo_1000 = np.array(geo_30)
            geo_1000[1] = geo_1000[1]*108
            geo_1000[5] = geo_1000[5]*108
            geo_1000 = tuple(geo_1000)
            
            # Create MEM up file
            dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
            dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)   
         
        BIO_map_Comp = os.path.join(composite_folder, date_comp, "biomass_prod_kg-ha_%s.tif" %(date_comp)) 
        dest_comp_up = RC.reproject_dataset_example(BIO_map_Comp, dest_up, 4) 
        
        # Downscale back to 30m
        dest_comp_down = RC.reproject_dataset_example(dest_comp_up, dest_LS_30, 2) 
        Array_comp = dest_comp_down.GetRasterBand(1).ReadAsArray()
                                  
        Temp_array = relMap_LS.clip(0.0001, 6) * Array_comp
        
        # Scale Temp_array
        dest_mem_temp = DC.Save_as_MEM(Temp_array, geo_30, 4326)
        
        # Open input composite and get geotransform
        dest_comp = gdal.Open(BIO_map_Comp)
        geo_comp = dest_comp.GetGeoTransform()
        proj_comp = dest_comp.GetProjection()
        Array_comp_MODIS = dest_comp.GetRasterBand(1).ReadAsArray()
        Array_comp_MODIS_BOX = Create_Relative_Map(Array_comp_MODIS, 3, Include_NAN = 1, type_return = "mean")
        
        # Upscale LS comp
        dest_LScomp_up = RC.reproject_dataset_example(dest_mem_temp, BIO_map_Comp, 4)
        Array_LScomp_up = dest_LScomp_up.GetRasterBand(1).ReadAsArray()
        Array_LScomp_up_BOX = Create_Relative_Map(Array_LScomp_up, 3, Include_NAN = 1, type_return = "mean")
        
        # Calculate Fraction Correction
        Corr_Frac = Array_comp_MODIS_BOX/Array_LScomp_up_BOX
        Corr_Frac = Corr_Frac.clip(0.0001, 6)
        dest_corr_up = DC.Save_as_MEM(Corr_Frac, geo_comp, proj_comp)
        
        # downscale corr factor
        dest_corr_down = RC.reproject_dataset_example(dest_corr_up, dest_mem_temp, 2)
        Corr_down = dest_corr_down.GetRasterBand(1).ReadAsArray()
        
        # apply corr
        End_array = Corr_down * Temp_array
        End_array[End_array>300] = 300
        
        DC.Save_as_tiff(filename_end_out, End_array, geo_30, 4326)
    
    return()
                            
            
        
def main_sm(Date, LS_folder, composite_folder, output_folder):

    date_comp = Date.strftime("%Y%m%d")
    filename_end_out = os.path.join(output_folder, "se_root_%s.tif" %(date_comp))   
 
    if not os.path.exists(filename_end_out):
        
        # Find all LS dates
        LS_LVL3_maps = os.listdir(LS_folder)
        DATES_LS_LVL3 = np.array([datetime.datetime.strptime(k, "%Y%m%d").toordinal() for k in LS_LVL3_maps if len(k) == 8])
        
        # Find amount of days between
        Diff_days = np.abs(DATES_LS_LVL3 - Date.toordinal())
        Diff_days_sorted = np.sort(Diff_days)
        
        # Create SE_rel map
        output_folder_SErel = os.path.join(LS_folder, "SErel")
        if not os.path.exists(output_folder_SErel):
            os.makedirs(output_folder_SErel)
            
        # 
        Continue = 1
        
        for Diff_day_sorted in [k for k in Diff_days_sorted if k <60]:
            
            if Continue == 1:
                
                take_days = np.argwhere(Diff_days == Diff_day_sorted)
                
                for take_day in take_days:
                
                    DATE_LS_LVL3 = int(DATES_LS_LVL3[take_day])
                    
                    date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")
                
                    print("Add %s in relative map" %date_LS)
                    
                    Output_filename_SErel = os.path.join(output_folder_SErel, "serel_%s.tif" %date_LS)
                    
                    if not os.path.exists(Output_filename_SErel):
                    
                        SE_map_LS = os.path.join(LS_folder, date_LS, "se_root_%s.tif" %(date_LS))
                        
                        dest_LS_30 = gdal.Open(SE_map_LS)
                        geo_30 = dest_LS_30.GetGeoTransform()
                        proj = dest_LS_30.GetProjection()
                        size_x = dest_LS_30.RasterXSize
                        size_y = dest_LS_30.RasterYSize
                        Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
                        
                        geo_1000 = np.array(geo_30)
                        geo_1000[1] = geo_1000[1]*108
                        geo_1000[5] = geo_1000[5]*108
                        geo_1000 = tuple(geo_1000)
                        
                        # Create MEM up file
                        dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
                        dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)
                        
                        # upscale LS
                        dest_LS_up = RC.reproject_dataset_example(dest_LS_30, dest_up, 4)
            
                        # Get relative array
                        relMap_LS_one = Create_Relative_Map(Array_LS, 108, Include_NAN = 0)  
                        
                        DC.Save_as_tiff(Output_filename_SErel, relMap_LS_one, geo_30, proj)
                    
                    else:
                        dest_serel = gdal.Open(Output_filename_SErel)
                        relMap_LS_one = dest_serel.GetRasterBand(1).ReadAsArray()
                    
                    if "relMap_LS" not in locals():
                        relMap_LS = np.ones(relMap_LS_one.shape) * np.nan
                        
                    # Fit in endmap
                    relMap_LS = np.where(np.isnan(relMap_LS), relMap_LS_one, relMap_LS)
                
                    Amount_Of_Pixels2Go = np.count_nonzero(np.isnan(relMap_LS))
                
                if Amount_Of_Pixels2Go == 0:
                    Continue = 0
                else:
                    print("Still %s pixels to go" %Amount_Of_Pixels2Go)
         
        if "dest_up" not in locals():
            
            DATE_LS_LVL3 = int(DATES_LS_LVL3[0])   
            date_LS = datetime.datetime.fromordinal(DATE_LS_LVL3).strftime("%Y%m%d")           
            SE_map_LS = os.path.join(LS_folder, date_LS, "se_root_%s.tif" %(date_LS))
                        
            dest_LS_30 = gdal.Open(SE_map_LS)
            geo_30 = dest_LS_30.GetGeoTransform()
            proj = dest_LS_30.GetProjection()
            size_x = dest_LS_30.RasterXSize
            size_y = dest_LS_30.RasterYSize
            Array_LS = dest_LS_30.GetRasterBand(1).ReadAsArray()
            
            geo_1000 = np.array(geo_30)
            geo_1000[1] = geo_1000[1]*108
            geo_1000[5] = geo_1000[5]*108
            geo_1000 = tuple(geo_1000)
            
            # Create MEM up file
            dummy_array = np.ones([int(np.floor(size_y/108)),int(np.floor(size_x/108))]) * np.nan
            dest_up = DC.Save_as_MEM(dummy_array, geo_1000, 4326)   
         
        SE_map_Comp = os.path.join(composite_folder, date_comp, "se_root_%s.tif" %(date_comp)) 
        dest_comp_up = RC.reproject_dataset_example(SE_map_Comp, dest_up, 4) 
        
        # Downscale back to 30m
        dest_comp_down = RC.reproject_dataset_example(dest_comp_up, dest_LS_30, 2) 
        Array_comp = dest_comp_down.GetRasterBand(1).ReadAsArray()
                                  
        Temp_array = relMap_LS.clip(0.0001, 6) * Array_comp
        
        # Scale Temp_array
        dest_mem_temp = DC.Save_as_MEM(Temp_array, geo_30, 4326)
        
        # Open input composite and get geotransform
        dest_comp = gdal.Open(SE_map_Comp)
        geo_comp = dest_comp.GetGeoTransform()
        proj_comp = dest_comp.GetProjection()
        Array_comp_MODIS = dest_comp.GetRasterBand(1).ReadAsArray()
        Array_comp_MODIS_BOX = Create_Relative_Map(Array_comp_MODIS, 3, Include_NAN = 1, type_return = "mean")
        
        # Upscale LS comp
        dest_LScomp_up = RC.reproject_dataset_example(dest_mem_temp, SE_map_Comp, 4)
        Array_LScomp_up = dest_LScomp_up.GetRasterBand(1).ReadAsArray()
        Array_LScomp_up_BOX = Create_Relative_Map(Array_LScomp_up, 3, Include_NAN = 1, type_return = "mean")
        
        # Calculate Fraction Correction
        Corr_Frac = Array_comp_MODIS_BOX/Array_LScomp_up_BOX
        Corr_Frac = Corr_Frac.clip(0.0001, 6)
        dest_corr_up = DC.Save_as_MEM(Corr_Frac, geo_comp, proj_comp)
        
        # downscale corr factor
        dest_corr_down = RC.reproject_dataset_example(dest_corr_up, dest_mem_temp, 2)
        Corr_down = dest_corr_down.GetRasterBand(1).ReadAsArray()
        
        # apply corr
        End_array = Corr_down * Temp_array
        End_array=End_array.clip(0,1)
        
        DC.Save_as_tiff(filename_end_out, End_array, geo_30, 4326)
    
    return()
                            
            
            
                   
        
    
   





'''
dest_out = RC.reproject_dataset_example(r"D:\Project_FAO\Test_LVL3\ET_LSCOMP_30m.tif", ET_map_Comp, 4)
array_out = dest_out.GetRasterBand(1).ReadAsArray()
dest_ex = gdal.Open(ET_map_Comp)
geo_ex = dest_ex.GetGeoTransform()
DC.Save_as_tiff(r"D:\Project_FAO\Test_LVL3\ET_LSCOMP_250m.tif", array_out, geo_ex, 4326)

Array_Comp = dest_comp_up.GetRasterBand(1).ReadAsArray()
Array_LS = dest_LS_up.GetRasterBand(1).ReadAsArray()

dest_LSCOMP_up = RC.reproject_dataset_example(r"D:\Project_FAO\Test_LVL3\ET_LSCOMP_30m.tif", dest_comp_up, 4)
Array_LSCOMP = dest_LSCOMP_up.GetRasterBand(1).ReadAsArray()
DC.Save_as_tiff(r"D:\Project_FAO\Test_LVL3\ET_Comp_1000m.tif", Array_Comp, geo_1000, 4326)
DC.Save_as_tiff(r"D:\Project_FAO\Test_LVL3\ET_LS_1000m.tif", Array_LS, geo_1000, 4326)
DC.Save_as_tiff(r"D:\Project_FAO\Test_LVL3\ET_LSCOMP_1000m.tif", Array_LSCOMP, geo_1000, 4326)


print("ET Average COMPOSITE = %s" %np.nanmean(Array_Comp))
print("ET Average LS = %s"%np.nanmean(Array_LS))
print("ET Average LS_COMPOSITE = %s"%np.nanmean(Array_LSCOMP))

'''










def Create_Relative_Map(Data_In, Box, Include_NAN = 1, type_return = "relative"):
    
    # Define required parameters
    Data_In[Data_In==-9999] = np.nan
    Size_box = int((Box-1)/2)
    Data_Out=np.zeros([len(Data_In),len(Data_In[1])]) 
    Data_amount = np.zeros([len(Data_In),len(Data_In[1])]) 

    # create 3D box of the moving window
    for ypixel in list(range(-Size_box,Size_box + 1)):
        for xpixel in list(range(-Size_box,Size_box + 1)):
           
            Xmin = np.maximum(0, xpixel)
            Xmax = np.minimum(len(Data_In[1]) + xpixel, len(Data_In[1]))
            Ymin = np.maximum(0, ypixel)
            Ymax = np.minimum(len(Data_In) + ypixel, len(Data_In))
            
            Xmin2 = len(Data_In[1]) - Xmax
            Xmax2 = len(Data_In[1]) - Xmin
            Ymin2 = len(Data_In) - Ymax
            Ymax2 = len(Data_In) - Ymin
            
            if Include_NAN == 1:
                Data_in_one = Data_In[Ymin:Ymax, Xmin:Xmax]
                Data_Out[Ymin2:Ymax2, Xmin2:Xmax2] += np.where(np.isnan(Data_in_one), 0, Data_in_one)
                
                Data_amount_one = np.where(np.isnan(Data_in_one), 0, 1)
                Data_amount[Ymin2:Ymax2, Xmin2:Xmax2]  += Data_amount_one
            else:
                Data_Out[Ymin2:Ymax2, Xmin2:Xmax2] += Data_In[Ymin:Ymax, Xmin:Xmax]
                Data_amount[Ymin2:Ymax2, Xmin2:Xmax2] += 1 
                     
    Data_mean = Data_Out/Data_amount

    if type_return == "relative":
        # Calculate the relative map
        End = Data_In/Data_mean
    else:
        End = Data_mean

    return(End)
