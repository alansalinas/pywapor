import os
import glob
import tqdm
import pyWAPOR
import pyWAPOR.tests.general as gen

init_file = os.path.abspath(pyWAPOR.__file__)
init_path = os.path.split(init_file)[0]
all_files = glob.glob(os.path.join(init_path, "**/*.py"), recursive=True)

iterations = tqdm.tqdm(all_files, 
                       desc = "Evaluating pylint-scores", 
                       unit = " files", 
                       position = 0, 
                       leave = True)

results = [gen.evaluate_pylint(file) for file in iterations]
output_fn = "pylint_scores_{0}_{1}.csv".format(pyWAPOR.__name__, 
                                               pyWAPOR.__version__)
output_path = os.path.join(init_path, "tests", "test_results", output_fn)
gen.create_csv(output_path, ['filename', 'pylint_score'], results)