import pyWAPOR
import pandas as pd
import os
import shutil
import gdal
import numpy as np

def ETLook_main(asserts = False):

    data_folder = os.path.join(pyWAPOR.__path__[0], "tests", "test_data")
    input_folder = os.path.join(data_folder, "input")
    output_folder = os.path.join(data_folder, "output")
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    date_string = "20190704"
    date = pd.Timestamp(date_string, freq="d")

    pyWAPOR.ETLook.ETLook_code.main(input_folder, output_folder, date)

    if asserts:

        vars = ["biomass_prod_kg-ha", "e_24_mm", "et_24_mm",
        "et_ref_24_mm", "int_mm", "se_root", "stress_vpd", "t_24_mm"]

        for var in vars:

            print(var)

            test_file = os.path.join(output_folder,
                                    date_string,
                                    "{0}_{1}.tif".format(var, date_string))

            output_ref_folder = os.path.join(data_folder, "output_ref")
            ref_file = os.path.join(output_ref_folder,
                                    date_string,
                                    "{0}_{1}.tif".format(var, date_string))

            test_ds = gdal.Open(test_file)
            ref_ds = gdal.Open(ref_file)

            test_array = test_ds.GetRasterBand(1).ReadAsArray()
            ref_array = ref_ds.GetRasterBand(1).ReadAsArray()

            test_nans = np.sum(np.isnan(test_array))
            ref_nans = np.sum(np.isnan(ref_array))

            diff_sum = np.nansum(np.abs(test_array - ref_array))

            assert test_nans == ref_nans, "Output has different amount of NaNs: {0} != {1}".format(test_nans, ref_nans)
            assert np.isclose(diff_sum, 0.0), "Output has different values: {0} != {1}".format(diff_sum, 0.0)

ETLook_main(asserts = True)