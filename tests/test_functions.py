#%%
import pywapor
import pandas as pd
import os
import shutil
import numpy as np
import pathlib
import pywapor.general.processing_functions as PF
import pywapor.general.variables as vars
import datetime

def compare_frame_pywapor(verbose = True):

    print("")

    minimal_keys = ['r0', 'ndvi', 'lst', 'dtime', 'lat_deg', 'z', 'land_mask', 'lue_max', 'p_air_0_24', 'p_air_0_i', 'p_air_i', 'P_24', 'qv_24', 'qv_i', 't_air_24', 't_air_i', 't_amp_year', 'u_24', 'u_i', 'wv_i', 'trans_24', 'se_root', 'z_obst_max', 'rs_min', 'ra_24', 'z_oro', "lw_offset", "lw_slope", "rn_offset", "rn_slope", "t_opt", "vpd_slope", "r0_bare", "r0_full"]

    data_fh = os.path.join(pathlib.Path(pywapor.__path__[0]).parent, 
                                "tests", "test_data", "wapor_frame_ref", 
                                "frame_test_sample_points_20210701.csv")

    df = pd.read_csv(data_fh)
    date = datetime.date(2021,7,1)
    names = [("1", "Agri. #1 "), ("2",  "Water    "),
            ("3",  "Urban    "), ("4",  "Desert   "),
            ("5",  "Agri. #2 ")
    ]

    input_data = dict()
    for param in minimal_keys:
        value = df.loc[df["parameter"] == param, [x[0] for x in names]].values
        if value.size == 0:
            input_data[param] = None
        else:
            input_data[param] = value

    od = pywapor.et_look.main("", date, et_look_version = "v2", input_data = input_data)[0]

    errors = np.array([[np.nan, np.nan, np.nan, np.nan, np.nan]])
    checks = np.array([]).astype(np.bool)
    checked_vars = np.array([])

    if not verbose:
        print("pixel_name".ljust(20) + f" --> {'; '.join([x[1] for x in names])}")

    for param in od.keys():
        
        ref_values = df.loc[df["parameter"] == param, [x[0] for x in names]].values
        calc_values = od[param]

        if ref_values.size == 0:
            continue

        error = np.abs(ref_values - calc_values)
        check = np.all(np.isclose(ref_values, calc_values, rtol = 1e-1))

        errors = np.concatenate((errors, error))
        checks = np.append(checks, check)
        checked_vars = np.append(checked_vars, param)

        if not verbose:
            if check:
                print(f"{param}".ljust(20) + f" --> {check}")
            else:
                print(f"{param}".ljust(20) + f" --> {'; '.join([f'{x:.7f}' for x in error[0]])}")

    errors = errors[1:,:]

    if verbose:

        check_e_24_mm = checks[checked_vars == "e_24_mm"][0]
        check_t_24_mm = checks[checked_vars == "t_24_mm"][0]
        
        if np.all([check_e_24_mm, check_t_24_mm]):
            print(f"E and T match between FRAME and pyWaPOR.")
        else:
            print(f"E and T do NOT match between FRAME and pyWaPOR.")

    print("")

    return zip(checked_vars, checks, errors)

def test_et_look_main(et_look_version, verbose = True):

    print("")

    data_folder = os.path.join(pathlib.Path(pywapor.__path__[0]).parent, 
                                "tests", "test_data")
    
    project_folder = os.path.join(data_folder, "input_ref")
    output_folder = os.path.join(project_folder, "out_level_1")
    
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    date_string = "20190704"
    date = pd.Timestamp(date_string, freq="d")

    pywapor.et_look.main(project_folder, date, et_look_version = et_look_version)

    nses = np.array([])
    checked_vars = np.array([])

    for var in vars.order:

        test_file = os.path.join(output_folder,
                                date_string,
                                "{0}_{1}.tif".format(var, date_string))

        output_ref_folder = os.path.join(data_folder, 
                                "output_ref_{0}".format(et_look_version))

        ref_file = os.path.join(output_ref_folder,
                                date_string,
                                "{0}_{1}.tif".format(var, date_string))

        if not np.all([os.path.isfile(test_file), os.path.isfile(ref_file)]):
            continue

        test_array = PF.open_as_array(test_file)
        ref_array = PF.open_as_array(ref_file)

        arrays = [test_array, ref_array]
        mask = np.any([np.isnan(array) for array in arrays], axis = 0)
        arrays = np.array([array[~mask] for array in arrays])

        nse = pywapor.post_et_look.calc_nash_sutcliffe(arrays)[0]

        nses = np.append(nses, nse)
        checked_vars = np.append(checked_vars, var)

        if not verbose:
            print(f"{var}".zfill(20).replace("0", " ") + f" --> {nse:.8f}")

    if verbose:
        if np.all(np.isclose(nses, np.ones_like(nses))):
            print(f"{et_look_version} working.")
        else:
            print(f"{et_look_version} NOT working.")

    print("")

    return zip(checked_vars, nses)

errors_frame_pywapor = compare_frame_pywapor(verbose = True)
nses_dev = test_et_look_main("dev")
nses_v2 = test_et_look_main("v2")
# %%
