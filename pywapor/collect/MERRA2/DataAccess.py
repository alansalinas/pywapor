
import os
import numpy as np
import pandas as pd
import datetime
import requests
from netCDF4 import Dataset
import tqdm
from pywapor.general.logger import log
import time

def DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, Period, username, password, Waitbar, data_type = ["mean"]):

    import pywapor.general.processing_functions as PF

    all_files = list()
    # if "mean" in data_type:
    #     all_files[Var] = list()
    # if "max" in data_type:
    #     all_files[f"{Var}-max"] = list()
    # if "min" in data_type:
    #     all_files[f"{Var}-min"] = list()

    # Add extra buffer to ensure good spatial interpolation
    buffer_pixels = 0
    lonlim = [lonlim[0] - 0.625 * buffer_pixels, lonlim[1] + 0.625 * buffer_pixels]
    latlim = [latlim[0] - 0.5 * buffer_pixels, latlim[1] + 0.5 * buffer_pixels]

    # Check the latitude and longitude and otherwise set lat or lon on greatest extent
    if latlim[0] < -90 or latlim[1] > 90:
        log.warning('Latitude above 90N or below 90S is not possible. Value set to maximum')
        latlim[0] = np.max(latlim[0], -90)
        latlim[1] = np.min(latlim[1], 90)
    if lonlim[0] < -180 or lonlim[1] > 180:
        log.warning('Longitude must be between 180E and 180W. Now value is set to maximum')
        lonlim[0] = np.max(lonlim[0], -180)
        lonlim[1] = np.min(lonlim[1], 180)  
    
    # Get information of the parameter    
    VarInfo = VariablesInfo(TimeStep)
    Parameter = VarInfo.names[Var]
    unit  = VarInfo.units[Var]
    p_ts = VariablesInfo.period_times

    # Create output folder
    output_folder = os.path.join(Dir, "MERRA2", Parameter, TimeStep) 
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    corrx = 0.625 * 0.5
    corry = 0.5 * 0.5

    # Define IDs
    IDx = [np.floor((lonlim[0] + corrx + 180)/0.625), np.ceil((lonlim[1] + corrx + 180)/0.625)]
    IDy = [np.floor((latlim[0] + corry + 90)/0.5), np.ceil((latlim[1] + corry + 90)/0.5)]
    
    # Create output geo transform
    Xstart = -180 + 0.625 * IDx[0] - corrx
    Ystart = -90 + 0.5 * IDy[1] - corry
    
    geo_out = tuple([Xstart, 0.625, 0, Ystart, 0, -0.5])
    proj = "WGS84"
    
    Dates = pd.date_range(Startdate, Enddate, freq = "D")

    if Waitbar:
        waitbar = tqdm.tqdm(desc= f"Tile: 0 / {len(Dates)}",
                            position = 0,
                            # total=total_size,
                            unit='Bytes',
                            unit_scale=True,)
    else:
        waitbar = None

    for Date in Dates:
        
        # Define the IDz
        if TimeStep == "hourly_MERRA2":
            if Period < 1 or Period > 24:
                raise ValueError(f"Invalid Period (={Period}).")

            date = datetime.datetime.combine(Date, 
                                    p_ts[Period])
            
            fn = f"{Var}_MERRA2_{unit}_inst_{date:%Y.%m.%d.%H.%M}.tif"
            output_name = os.path.join(output_folder, fn)

            output_folder_temp = os.path.join(Dir, "MERRA2", "Temp")
            if not os.path.exists(output_folder_temp):
                os.makedirs(output_folder_temp)
            year = Date.year
            month = Date.month
            day = Date.day
            output_name_min = output_folder
            output_name_max = output_folder
            if "min" in data_type:
                log.warning(f"data_type = {data_type} not applicable to hourly_MERRA2, ignoring.")
                data_type.remove("min")
            if "max" in data_type:
                log.warning(f"data_type = {data_type} not applicable to hourly_MERRA2, ignoring.")
                data_type.remove("max")
            
        if TimeStep == "daily_MERRA2":

            if "mean" in data_type:
                output_name = os.path.join(output_folder, "%s_MERRA2_%s_daily_%d.%02d.%02d.tif" %(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name = output_folder
            if "min" in data_type:
                output_name_min = os.path.join(output_folder, "%smin_MERRA2_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name_min = output_folder
            if "max" in data_type:
                output_name_max = os.path.join(output_folder, "%smax_MERRA2_%s_daily_%d.%02d.%02d.tif"%(Var, unit, Date.year, Date.month, Date.day))
            else:
                output_name_max = output_folder
                
            output_folder_temp = os.path.join(Dir, "MERRA2", "Temp")
            if not os.path.exists(output_folder_temp):
                os.makedirs(output_folder_temp)            
            year = Date.year
            month = Date.month
            day = Date.day
    
        if not (os.path.exists(output_name) and os.path.exists(output_name_min) and os.path.exists(output_name_max)):

            # TODO this shouldnt be hardcoded, cause it will keep giving problems going forward.
            if Date < datetime.datetime(1992,1,1):
                number = 1
            elif (Date >= datetime.datetime(1992,1,1) and Date < datetime.datetime(2001,1,1)):
                number = 2
            elif (Date >= datetime.datetime(2001,1,1) and Date < datetime.datetime(2011,1,1)):
                number = 3
            else:
                number = 4
            if Date.month == 9 and Date.year == 2020:
                number2 = 1
            elif Date.month == 6 and Date.year == 2021:
                number2 = 1
            elif Date.month == 7 and Date.year == 2021:
                number2 = 1
            elif Date.month == 8 and Date.year == 2021:
                number2 = 1
            elif Date.month == 9 and Date.year == 2021:
                number2 = 1
            else:
                number2 = 0
                            
            if Var == "swgnet":
                url_MERRA = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXRAD.5.12.4/%d/%02d/MERRA2_%s0%s.tavg1_2d_rad_Nx.%d%02d%02d.nc4" %(year, month, number, number2, year, month, day)
            else:    
                url_MERRA = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/%d/%02d/MERRA2_%s0%s.inst1_2d_asm_Nx.%d%02d%02d.nc4" %(year, month, number, number2, year, month, day)
                
            # Reset the begin parameters for downloading
            downloaded = 0
            N = 0

            # if not downloaded try to download file
            while downloaded == 0:

                # Define the output name that is downloaded
                file_name = os.path.join(output_folder_temp, url_MERRA.split("/")[-1])

                if os.path.exists(file_name):
                    done = os.stat(file_name).st_size > 1000
                    if not done:
                        os.remove(file_name)
                else:
                    done = False

                if not done:
                
                    # make contact with server
                    x = requests.get(url_MERRA, allow_redirects = False)

                    x.raise_for_status()

                    resp = requests.get(x.headers['location'], auth = (username, password), stream=True)
                    total_size = int(resp.headers.get('content-length', 0))
                    
                    try:
                        resp.raise_for_status()
                    except requests.exceptions.HTTPError as e:
                        if N < 5:

                            N += 1
                            log.warning(f"N = {N}, {url_MERRA}")
                            time.sleep(5)
                            continue
                        else:
                            raise e

                    if not isinstance(waitbar, type(None)):
                        waitbar.reset(total = total_size)
                        waitbar_i = int(waitbar.desc.split(" ")[1])
                        waitbar_desc = str(waitbar.desc)
                        waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /"))
                    with open(file_name, 'wb') as z:
                        for data in resp.iter_content(chunk_size=1024):
                            size = z.write(data)
                            if not isinstance(waitbar, type(None)):
                                waitbar.update(size)
                                
                    statinfo = os.stat(file_name)
                    # Say that download was succesfull
                    if int(statinfo.st_size) > 1000:
                            downloaded = 1  

                else:
                    downloaded = 1
                    if not isinstance(waitbar, type(None)):
                        waitbar_i = int(waitbar.desc.split(" ")[1])
                        waitbar_desc = str(waitbar.desc)
                        waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /")) 

            data_end, data_min, data_max = Get_NC_data_end(file_name,Var, TimeStep, Period, IDy, IDx, VarInfo)

            # Add the VarFactor
            if VarInfo.factors[Var] < 0:
                if "mean" in data_type:                       
                    data_end[data_end != -9999] = data_end[data_end != -9999] + VarInfo.factors[Var]
                    data_end[data_end < -9999] = -9999
                    data_end = np.flipud(data_end)
                if "min" in data_type:
                    data_min[data_min != -9999] = data_min[data_min != -9999] + VarInfo.factors[Var]
                    data_min[data_min < -9999] = -9999
                    data_min = np.flipud(data_min)                                 
                if "max" in data_type:
                    data_max[data_max != -9999] = data_max[data_max != -9999] + VarInfo.factors[Var]      
                    data_max[data_max < -9999] = -9999
                    data_max = np.flipud(data_max)                        
                                            
            else:
                if "mean" in data_type: 
                    data_end[data_end != -9999] = data_end[data_end != -9999] * VarInfo.factors[Var]
                    data_end[data_end < -9999] = -9999
                    data_end = np.flipud(data_end)
                if "min" in data_type:
                    data_min[data_min != -9999] = data_min[data_min != -9999] * VarInfo.factors[Var]
                    data_min[data_min < -9999] = -9999
                    data_min = np.flipud(data_min)                            
                if "max" in data_type:
                    data_max[data_max != -9999] = data_max[data_max != -9999] * VarInfo.factors[Var]                          
                    data_max[data_max < -9999] = -9999
                    data_max = np.flipud(data_max)                        
                
            # Save as tiff file
            if "mean" in data_type: 
                PF.Save_as_tiff(output_name, data_end, geo_out, proj)
                all_files.append(output_name)
            if "min" in data_type:
                PF.Save_as_tiff(output_name_min, data_min, geo_out, proj)
                all_files.append(output_name_min)
            if "max" in data_type:
                PF.Save_as_tiff(output_name_max, data_max, geo_out, proj)
                all_files.append(output_name_max)

        else:
            if os.path.isfile(output_name):
                all_files.append(output_name)
            if os.path.isfile(output_name_min):
                all_files.append(output_name_min)
            if os.path.isfile(output_name_max):
                all_files.append(output_name_max)

            if not isinstance(waitbar, type(None)):
                waitbar_i = int(waitbar.desc.split(" ")[1])
                waitbar_desc = str(waitbar.desc)
                waitbar.set_description_str(waitbar_desc.replace(f": {waitbar_i} /", f": {waitbar_i+1} /")) 

    return all_files  

def Get_NC_data_end(file_name,Var, TimeStep, Period, IDy, IDx, VarInfo):
                        
    dict_para = {'t2m': 'T2M',
         'u2m': 'U2M',
         'v2m': 'V2M',
         'q2m': 'QV2M',
         'tpw': 'TQV',
         'ps': 'PS',
         'slp': 'SLP',
         'swgnet': 'SWGDN'}  
      
    types  = VarInfo.types[Var]
    if TimeStep == "hourly_MERRA2":
        data_end = Dataset(file_name)["%s" %dict_para[Var]][int(Period-1),int(IDy[0]):int(IDy[1]),int(IDx[0]):int(IDx[1])]
        data_min = []
        data_max = []
        
    else:
        data = Dataset(file_name)["%s" %dict_para[Var]][:,int(IDy[0]):int(IDy[1]),int(IDx[0]):int(IDx[1])]
        if types == "state":
            data_end = np.nanmean(data, 0)
        else:
            data_end = np.nansum(data, 0)
            
        data[data==-9999] = np.nan
        data_min = np.nanmin(data, 0)
        data_max = np.nanmax(data, 0)     
        
    return(data_end, data_min, data_max)

class VariablesInfo:
    """
    This class contains the information about the GLDAS variables
    """
    names = {'t2m': 'Air_Temperature',
             'u2m': 'Eastward_Wind',
             'v2m': 'Northward_Wind',
             'q2m': 'Specific_Humidity',
             'tpw': 'Total_Precipitable_Water_Vapor',
             'ps': 'Surface_Pressure',
             'slp': 'Sea_Level_Pressure',
             'swgnet': 'Surface_Net_Downward_Shortwave_Flux'
             }
    
    descriptions = {'t2m': '2m Air Temperature',
             'u2m': '2m Eastward wind',
             'v2m': '2m Northward wind',
             'q2m': '2m Specific Humidity',
             'tpw': 'Total Precipitable Water Vapor',
             'ps': 'Surface Pressure',
             'slp': 'Sea Level Pressure',
             'swgnet': 'Surface Net Downward Shortwave Flux'
             }
    
    factors = {'t2m': 1,
             'u2m': 1,
             'v2m': 1,
             'q2m': 1,
             'tpw': 1,
             'ps': 0.001,
             'slp':  0.001,
             'swgnet':  1
             }
    
    types = {'t2m': 'state',
             'u2m': 'state',
             'v2m': 'state',
             'q2m': 'state',
             'tpw': 'state',
             'ps': 'state',
             'slp': 'state',
             'swgnet': 'state'
             }

    period_times = {1: datetime.time(0, 30),
                    2: datetime.time(1, 30),
                    3: datetime.time(2, 30),
                    4: datetime.time(3, 30),
                    5: datetime.time(4, 30),
                    6: datetime.time(5, 30),
                    7: datetime.time(6, 30),
                    8: datetime.time(7, 30),
                    9: datetime.time(8, 30),
                    10: datetime.time(9, 30),
                    11: datetime.time(10, 30),
                    12: datetime.time(11, 30),
                    13: datetime.time(12, 30),
                    14: datetime.time(13, 30),
                    15: datetime.time(14, 30),
                    16: datetime.time(15, 30),
                    17: datetime.time(16, 30),
                    18: datetime.time(17, 30),
                    19: datetime.time(18, 30),
                    20: datetime.time(19, 30),
                    21: datetime.time(20, 30),
                    22: datetime.time(21, 30),
                    23: datetime.time(22, 30),
                    24: datetime.time(23, 30)}

    def __init__(self, step):
            
        if step == 'hourly_MERRA2':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }    
            
        elif step == 'daily_MERRA2':
            self.units = {'t2m': 'K',
             'u2m': 'm-s-1',
             'v2m': 'm-s-1',
             'q2m': 'kg-kg-1',
             'tpw': 'mm',
             'ps': 'kpa',
             'slp': 'kpa',
             'swgnet': 'W-m-2'
             }

        else:
            raise KeyError("The input time step is not supported")

if __name__ == "__main__":

    import pywapor

    Dir = r"/Users/hmcoerver/pywapor_notebooks"
    Startdate = "2021-06-01"
    Enddate = "2021-06-02"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    username, password = pywapor.collect.get_pw_un.get("NASA")
    
    Period = 1
    
    data_type = ["mean", "min", "max"]

    print("daily")
    TimeStep = "daily_MERRA2"
    meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp', 'swgnet']
    # meteo_vars = ["swgnet"]
    for Var in meteo_vars:
        print(Var)
        all_files1 = DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, "", username, password, Waitbar = True, data_type = data_type)

    # print("hourly")
    # TimeStep = "hourly_MERRA2"
    # meteo_vars = ['t2m', 'u2m', 'v2m', 'q2m', 'tpw', 'ps', 'slp']
    # data_type = ["mean"]
    # for Var in meteo_vars:
    #     print(Var)
    #     all_files2 = DownloadData(Dir, Var, Startdate, Enddate, latlim, lonlim, TimeStep, Period, username, password, Waitbar = True, data_type = data_type)

