import pyWAPOR
import pandas as pd
import os

def profiler_ETLook(data_folder):

    input_folder = os.path.join(data_folder, "input")
    output_folder = os.path.join(data_folder, "output")
    date = pd.Timestamp('2019-07-04 00:00:00', freq="d")

    pyWAPOR.ETLook.ETLook_code.main(input_folder, output_folder, date)
