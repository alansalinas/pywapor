import cProfile
import os
import pyWAPOR
import pyWAPOR.tests.test_functions
import pyWAPOR.tests.general as gen

data_folder = os.path.join(pyWAPOR.__path__[0], r"tests/test_data")

function_path = ["pyWAPOR", "tests", "test_functions", "profiler_ETLook"]      
function_inputs = "('{0}')".format(data_folder)
command = ".".join(function_path) + function_inputs
stats_file = os.path.join(data_folder, "stats")

cProfile.run(command, stats_file)
gen.process_cProfile_stats(stats_file, version = pyWAPOR.__version__)