import glob
import os
import pywapor
import watertools
import datetime

# inputs
latlim = [26.06, 27.6]
lonlim = [78.55, 80.6]
output_folder = r"/Volumes/Data/RS_scenes/india_case/et_look_run"
Satellite_folder = r"/Volumes/Data/RS_scenes/india_case/C1_L1_output/145041"

tar_files = glob.glob(os.path.join(r"/Volumes/Data/RS_scenes/india_case/C1_L1", "*.tar.gz"))
dates = [datetime.datetime.strptime(k.split("_")[-4], "%Y%m%d") for k in tar_files]
Startdate = sorted(dates)[0].strftime("%Y-%m-%d")
Enddate = sorted(dates)[-1].strftime("%Y-%m-%d")

# Create 8 daily timesteps
dates_8d = watertools.Collect.MOD11.DataAccess.Make_TimeStamps(Startdate, Enddate)

# Collect data for the whole daily MODIS period
pywapor.pre_et_look.main(output_folder,
                        Startdate,
                        Enddate,
                        latlim,
                        lonlim,
                        LandCover = "GlobCover",
                        Short_Downwards_Radiation = "MERRA",
                        composite = True,
                        RAW_folder=os.path.join(output_folder,"RAW"))

# Run pyWAPOR model over daily MODIS data
for date in dates_8d:
    pywapor.et_look_code.main(os.path.join(output_folder, "ETLook_input_composite"),
                                    os.path.join(output_folder, "ETLook_output_composite"),
                                    date)

for Date in dates:

    Startdate = Date.strftime("%Y-%m-%d")
    Enddate = Date.strftime("%Y-%m-%d")

    # Collect date for the LS day
    pywapor.pre_et_look.main(os.path.join(output_folder),
                            Startdate,
                            Enddate,
                            latlim,
                            lonlim,
                            LandCover = "GlobCover",
                            Short_Downwards_Radiation = "MERRA",
                            Satellite_folder = Satellite_folder,
                            composite = False,
                            RAW_folder=os.path.join(output_folder,"RAW"))

    # Run pyWAPOR model for LS day
    pywapor.et_look_code.main(os.path.join(output_folder, "ETLook_input"),
                              os.path.join(output_folder, "ETLook_output_newconstants"),
                                    Date)

# Combine MODIS composite run and LS daily runs
for date in dates_8d:
    print(date)
    pywapor.calculate_composite_ls_modis.main_et(date,
                                                 os.path.join(output_folder, "ETLook_output"),
                                                 os.path.join(output_folder, "ETLook_output_composite"),
                                                 os.path.join(output_folder, "ETLook_output_composite_LS_MODIS"))