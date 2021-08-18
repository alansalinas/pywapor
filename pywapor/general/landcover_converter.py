# -*- coding: utf-8 -*-
"""
WaterSat
author: Tim Martijn Hessels
Created on Tue Feb 26 08:17:10 2019
"""

def Globcover_LM(version = '1.0'):
    ETlook_LM = {
    11: 1,   #Post-flooding or irrigated croplands (or aquatic)
    14:	1,   #Rainfed croplands
    20:	1,   #Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
    30:	1,   #Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
    40:	1,   #Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
    50:	1,   #Closed (>40%) broadleaved deciduous forest (>5m)
    60:	1,   #Open (15-40%) broadleaved deciduous forest/woodland (>5m)
    70:	1,   #Closed (>40%) needleleaved evergreen forest (>5m)
    90:	1,   #Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
    100: 1,   	#Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
    110: 1,   	#Mosaic forest or shrubland (50-70%) / grassland (20-50%)
    120: 1,   	#Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
    130: 1,  	#Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
    140: 1,  	#Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
    150: 1,  	#Sparse (<15%) vegetation
    160: 1,  	#Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
    170: 1,  	#Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
    180: 1,  	#Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
    190: 3,  	#Artificial surfaces and associated areas (Urban areas >50%)
    200: 1,  	#Bare areas
    210: 2,  	#Water bodies
    220: 1,  	#Permanent snow and ice
    230: 0  	#No data (burnt areas, clouds,…)
    }

    Classes_LM =dict()
    Classes_LM['1.0'] = ETlook_LM

    return Classes_LM[version]


def Globcover_MaxObs(version = '1.0'):
    ETlook_Classes = {
    11: 4.0,   #Post-flooding or irrigated croplands (or aquatic)
    14:	4.0,   #Rainfed croplands
    20:	2.0,   #Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
    30:	3.5,   #Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
    40:	0.1,   #Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
    50:	0.6,   #Closed (>40%) broadleaved deciduous forest (>5m)
    60:	1.2,   #Open (15-40%) broadleaved deciduous forest/woodland (>5m)
    70:	2.0,   #Closed (>40%) needleleaved evergreen forest (>5m)
    90:	5.0,   #Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
    100: 8.0,   	#Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
    110: 2.0,   	#Mosaic forest or shrubland (50-70%) / grassland (20-50%)
    120: 8.0,   	#Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
    130: 4.0,  	#Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
    140: 2.0,  	#Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
    150: 1.0,  	#Sparse (<15%) vegetation
    160: 0.3,  	#Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
    170: 6.0,  	#Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
    180: 10,  	#Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
    190: 10,  	#Artificial surfaces and associated areas (Urban areas >50%)
    200: 0.1,  	#Bare areas
    210: 0.1,  	#Water bodies
    220: 0.1,  	#Permanent snow and ice
    230: 0  	#No data (burnt areas, clouds,…)
    }

    Classes_MaxObs =dict()
    Classes_MaxObs['1.0'] = ETlook_Classes

    return Classes_MaxObs[version]


def Globcover_Bulk(version = '1.0'):
    ETlook_Classes = {
    11: 200,   #Post-flooding or irrigated croplands (or aquatic)
    14:	200,   #Rainfed croplands
    20:	150,   #Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
    30:	150,   #Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
    40:	100,   #Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
    50:	120,   #Closed (>40%) broadleaved deciduous forest (>5m)
    60:	100,   #Open (15-40%) broadleaved deciduous forest/woodland (>5m)
    70:	150,   #Closed (>40%) needleleaved evergreen forest (>5m)
    90:	180,   #Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
    100: 175,   	#Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
    110: 150,   	#Mosaic forest or shrubland (50-70%) / grassland (20-50%)
    120: 350,   	#Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
    130: 175,  	#Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
    140: 250,  	#Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
    150: 150,  	#Sparse (<15%) vegetation
    160: 250,  	#Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
    170: 200,  	#Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
    180: 300,  	#Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
    190: 100,  	#Artificial surfaces and associated areas (Urban areas >50%)
    200: 100,  	#Bare areas
    210: 100,  	#Water bodies
    220: 100,  	#Permanent snow and ice
    230: 0  	#No data (burnt areas, clouds,…)
    }

    Classes_Bulk =dict()
    Classes_Bulk['1.0'] = ETlook_Classes

    return Classes_Bulk[version]

def Globcover_LUEmax(version = '1.0'):
    ETlook_Classes = {
    11: 3.04,   #Post-flooding or irrigated croplands (or aquatic)
    14:	3.04,   #Rainfed croplands
    20:	3.04,   #Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
    30:	3.04,   #Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
    40:	2.45,   #Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
    50:	2.66,   #Closed (>40%) broadleaved deciduous forest (>5m)
    60:	2.66,   #Open (15-40%) broadleaved deciduous forest/woodland (>5m)
    70:	2.5,   #Closed (>40%) needleleaved evergreen forest (>5m)
    90:	2.5,   #Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
    100: 2.54,   	#Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
    110: 2.54,   	#Mosaic forest or shrubland (50-70%) / grassland (20-50%)
    120: 2.54,   	#Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
    130: 2.54,  	#Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
    140: 2.3,  	#Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
    150: 1.57,  	#Sparse (<15%) vegetation
    160: 2.66,  	#Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
    170: 1.57,  	#Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
    180: 1.57,  	#Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
    190: 1.57,  	#Artificial surfaces and associated areas (Urban areas >50%)
    200: 1.57,  	#Bare areas
    210: 0,  	#Water bodies
    220: 0,  	#Permanent snow and ice
    230: 0  	#No data (burnt areas, clouds,…)
    }

    Classes_LUEmax =dict()
    Classes_LUEmax['1.0'] = ETlook_Classes

    return Classes_LUEmax[version]



def WAPOR_MaxObs(version = '1.0'):
    ETlook_Classes = {
    33: 8.0,   #Natural Vegetation
    41:	2.0,   #Cropland, Rainfed or Fallow
    42:	4.0,   #Cropland, Irrigated or Under Water Management
    50:	10,   #Artificial Surfaces
    60:	0.1,   #Bare Soil Or Sparse Vegetation
    80:	0.1   #Water
    }

    Classes_MaxObs =dict()
    Classes_MaxObs['1.0'] = ETlook_Classes

    return Classes_MaxObs[version]

def WAPOR_Bulk(version = '1.0'):
    ETlook_Classes = {
    33: 175,   #Natural Vegetation
    41:	150,   #Cropland, Rainfed or Fallow
    42:	200,   #Cropland, Irrigated or Under Water Management
    50:	400,   #Artificial Surfaces
    60:	100,   #Bare Soil Or Sparse Vegetation
    80:	100   #Water
    }

    Classes_Bulk =dict()
    Classes_Bulk['1.0'] = ETlook_Classes

    return Classes_Bulk[version]


def WAPOR_LM(version = '1.0'):
    ETlook_Classes = {
    33: 1,   #Natural Vegetation
    41:	1,   #Cropland, Rainfed or Fallow
    42:	1,   #Cropland, Irrigated or Under Water Management
    50:	3,   #Artificial Surfaces
    60:	1,   #Bare Soil Or Sparse Vegetation
    80:	2   #Water
    }

    Classes_Bulk =dict()
    Classes_Bulk['1.0'] = ETlook_Classes

    return Classes_Bulk[version]


def WAPOR_LUEmax(version = '1.0'):
    ETlook_Classes = {
    33: 2.66,   #Natural Vegetation
    41:	3.04,   #Cropland, Rainfed or Fallow
    42:	3.04,   #Cropland, Irrigated or Under Water Management
    50:	1.57,   #Artificial Surfaces
    60:	2.66,   #Bare Soil Or Sparse Vegetation
    80:	0   #Water
    }

    Classes_LUEmax =dict()
    Classes_LUEmax['1.0'] = ETlook_Classes

    return Classes_LUEmax[version]

def get_lulc_values():
    lulc_values = {"GLOBCOVER": {"LandMask": Globcover_LM(),
                                "Bulk_Stomatal_resistance": Globcover_Bulk(),
                                "Maximum_Obstacle_Height": Globcover_MaxObs(),
                                "LUEmax": Globcover_LUEmax()},
                    "WAPOR":     {"LandMask": WAPOR_LM(),
                                "Bulk_Stomatal_resistance": WAPOR_Bulk(),
                                "Maximum_Obstacle_Height": WAPOR_MaxObs(),
                                "LUEmax": WAPOR_LUEmax()}
                            }

    return lulc_values