#%%
import pywapor
import pandas as pd
import os
import shutil
import numpy as np
import pathlib
import datetime
import pywapor.general.processing_functions as PF

def ETLook_main(asserts = False, ETLook_version = "dev"):

    data_folder = os.path.join(pathlib.Path(pywapor.__path__[0]).parent, 
                                "tests", "test_data")
    input_folder = os.path.join(data_folder, "input_ref")
    output_folder = os.path.join(data_folder, "output")
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    date_string = "20190704"
    date = pd.Timestamp(date_string, freq="d")

    pywapor.et_look_code.main(input_folder, output_folder, date, ETLook_version = ETLook_version)

    if asserts:

        vars = ["e_24_mm", "et_24_mm",
        "et_ref_24_mm", "int_mm", "se_root", "stress_vpd", "t_24_mm"]

        if ETLook_version == "dev":
            vars.append("biomass_prod_kg-ha")
  
        for var in vars:

            test_file = os.path.join(output_folder,
                                    date_string,
                                    "{0}_{1}.tif".format(var, date_string))

            output_ref_folder = os.path.join(data_folder, 
                                    "output_ref_{0}".format(ETLook_version))

            ref_file = os.path.join(output_ref_folder,
                                    date_string,
                                    "{0}_{1}.tif".format(var, date_string))

            test_array = PF.open_as_array(test_file)
            ref_array = PF.open_as_array(ref_file)

            test_nans = np.sum(np.isnan(test_array))
            ref_nans = np.sum(np.isnan(ref_array))

            diff_sum = np.nansum(np.abs(test_array - ref_array))

            assert test_nans == ref_nans, "Output has different amount of NaNs: {0} != {1}".format(test_nans, ref_nans)
            assert np.isclose(diff_sum, 0.0), "Output has different values: {0} != {1}".format(diff_sum, 0.0)

            print("CORRECT: {0}".format(var))

ETLook_main(asserts = True, ETLook_version = "v2")

# %%

def test_collect():

    data_folder = os.path.join(pathlib.Path(pywapor.__path__[0]).parent, 
                                "tests", "test_data")
    folders_input_RAW = os.path.join(data_folder, "RAW")
    if os.path.exists(folders_input_RAW):
        shutil.rmtree(folders_input_RAW)

    date = datetime.datetime(2019, 7, 6)

    startdate = date.strftime("%Y-%m-%d")
    enddate = date.strftime("%Y-%m-%d")
    latlim = [29.0, 29.6]
    lonlim = [30.3, 31.1]

    # un = passwords.passes["NASA"][0]
    # pw = passwords.passes["NASA"][1]
    # token = passwords.passes["WAPOR"][1]

    # un_vito = passwords.passes["VITO"][0]
    # pw_vito = passwords.passes["VITO"][1]

    # pywapor.collect.PROBAV.PROBAV_S5(folders_input_RAW, startdate, enddate, 
    # latlim, lonlim, un_vito, pw_vito, buffer_dates = False)

    # date = datetime.datetime(2002, 7, 4)
    # startdate = date.strftime("%Y-%m-%d")
    # enddate = date.strftime("%Y-%m-%d")    
    # pywapor.collect.MYD11.LST(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw)

    # date = datetime.datetime(2000, 2, 24)
    # startdate = date.strftime("%Y-%m-%d")
    # enddate = date.strftime("%Y-%m-%d")    
    # pywapor.Collect.MOD11.LST(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw, remove_hdf = 0)

    # date = datetime.datetime(2017, 12, 1) # first date
    # startdate = date.strftime("%Y-%m-%d")
    # enddate = date.strftime("%Y-%m-%d")     
    # pywapor.Collect.GEOS.daily(folders_input_RAW, ['t2m'],startdate, enddate, latlim, lonlim)
    # pywapor.Collect.GEOS.three_hourly(folders_input_RAW, ['u2m'], startdate, enddate, latlim, lonlim, [1])

    # date = datetime.datetime(1980, 1, 1)
    # startdate = date.strftime("%Y-%m-%d")
    # enddate = date.strftime("%Y-%m-%d") 
    # pywapor.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, un, pw)
    # pywapor.Collect.MERRA.hourly_MERRA2(folders_input_RAW, ['u2m'], startdate, enddate, latlim, lonlim, [1], un, pw)

    # pywapor.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, un, pw, data_type = ["mean", "min", "max"])
    # pywapor.Collect.MERRA.hourly_MERRA2(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, un, pw, [1])

    # date = datetime.datetime(1981, 1, 1)
    # startdate = date.strftime("%Y-%m-%d")
    # enddate = date.strftime("%Y-%m-%d") 
    # pywapor.Collect.CHIRPS.daily(folders_input_RAW, startdate, enddate, latlim, lonlim)

    # date = datetime.datetime(2009, 1, 1)
    # pywapor.Collect.WAPOR.Get_Layer(folders_input_RAW, date.strftime("%Y-01-01"), date.strftime("%Y-12-31"), latlim, lonlim, 'L1_LCC_A', token)

    # pywapor.collect.MOD13.NDVI(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw, remove_hdf = 0)
    # pywapor.Collect.MYD13.NDVI(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw)
    
    # pywapor.Collect.MCD43.ALBEDO(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw)
    # pywapor.Collect.MOD11.LST(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw)
    # pywapor.Collect.MYD11.LST(folders_input_RAW, startdate, enddate, latlim, lonlim, un, pw)

    # pywapor.Collect.CHIRPS.daily(folders_input_RAW, startdate, enddate, latlim, lonlim)
    # pywapor.Collect.SRTM.DEM(folders_input_RAW, latlim, lonlim)

    # pywapor.Collect.Globcover.Landuse(folders_input_RAW, latlim, lonlim)
    # pywapor.Collect.WAPOR.Get_Layer(folders_input_RAW, date.strftime("%Y-01-01"), date.strftime("%Y-12-31"), latlim, lonlim, 'L1_LCC_A', token)

    # pywapor.Collect.MERRA.yearly_T_Amplitude(folders_input_RAW, [date.year],latlim, lonlim)
    # pywapor.Collect.MERRA.daily(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim)
    # pywapor.Collect.MERRA.three_hourly(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, [1])
    # pywapor.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, un, pw)
    # pywapor.Collect.MERRA.hourly_MERRA2(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, [1], un, pw)
    # pywapor.Collect.MERRA.daily(folders_input_RAW, ['swgnet'], startdate, enddate, latlim, lonlim)
    # pywapor.Collect.MERRA.daily_MERRA2(folders_input_RAW, ['swgnet'],startdate, enddate, latlim, lonlim, un, pw)

    # pywapor.Collect.GEOS.daily(folders_input_RAW, ['t2m'],startdate, enddate, latlim, lonlim)
    # pywapor.Collect.GEOS.three_hourly(folders_input_RAW, ['t2m'], startdate, enddate, latlim, lonlim, [1])
    
    # pywapor.Collect.MSGCPP.SDS(folders_input_RAW, StartTime, EndTime, latlim, lonlim)

# test_collect()