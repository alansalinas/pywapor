import csv
import os
import pstats
from pylint.reporters import CollectingReporter
from pylint.lint import Run

def create_csv(output_path, header_row, data_rows):
    with open(output_path, 'w') as csvfile: 
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerow(header_row) 
        csvwriter.writerows(data_rows)

def split_path(path, splitter):
    path = os.path.normpath(path)
    parts = path.split(os.sep)
    idx = parts.index(splitter) + 1
    fn = os.path.join(*parts[idx:])
    return fn

def parse_info(x):
    if x[-1:] == "\n":
        path, rest = x.split(":")
        if "pyWAPOR" in path:
            path = split_path(path, "pyWAPOR")
        line_no, func_name = rest.split("(")
        func_name = func_name[:-2]
        return [path, line_no, func_name]
    else:
        return x

def flatten(lines):
    flat_lines = list()
    for row in lines:
        flat_sublines = list()
        for item in row:
            if isinstance(item, list):
                for subitem in item:
                    flat_sublines.append(subitem)
            else:
                flat_sublines.append(item)
        flat_lines.append(flat_sublines)
    return flat_lines

def evaluate_pylint(file):
    fn = split_path(file, "pyWAPOR")
    rep = CollectingReporter()
    files = [file, "-sn"]
    Run(files, reporter = rep, do_exit = False)
    fn_score = rep.linter.stats["global_note"]
    rep = None
    return [fn, "{0:.2f}".format(fn_score)]

def process_cProfile_stats(stats_file, version = "NaN"):
    stats_txt_file = stats_file + ".txt"
    
    with open(stats_txt_file, 'w') as stream:
        stats = pstats.Stats(stats_file, stream=stream)
        stats.sort_stats(pstats.SortKey.CUMULATIVE).print_stats("pyWAPOR")

    with open(stats_txt_file) as f:
        lines = f.readlines()

    os.remove(stats_txt_file)
    os.remove(stats_file)

    header_idx = [i for i, line in enumerate(lines) if "ncalls" in  line][0]
    clean_lines = [[parse_info(x) for x in line.split(" ") if x != ""] for line in lines[header_idx:-2]]
    flat_lines = flatten(clean_lines)

    create_csv("test_results/profiler_ETLook_{0}.csv".format(version), flat_lines[0], flat_lines[1:])