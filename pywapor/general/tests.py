import requests
import pywapor
import os
import time
import numpy as np

def check_source_selection(source_selection, startdate, enddate):

    valid_sources, valid_dates = pywapor.general.variables.get_source_validations()

    temporal_sources = ["METEO", "NDVI", "ALBEDO", "LST", 
                    "PRECIPITATION", "TRANS"]

    check_keys = np.all([key in valid_sources.keys() for key in source_selection.keys()])
    assert check_keys, "invalid key in source_selection"

    assert len(source_selection["DEM"]) == 1, "only one DEM source can be selected"
    assert len(source_selection["METEO"]) == 1, "only one METEO source can be selected"
    assert len(source_selection["PRECIPITATION"]) == 1, "only one PRECIPITATION source can be selected"
    assert len(source_selection["LULC"]) == 1, "only one LULC source can be selected"
    assert len(source_selection["TRANS"]) == 1, "only one TRANS source can be selected"

    results = dict()
    all_results = list()

    for var, sources in source_selection.items():

        check1 = [source in valid_sources[var] for source in sources]
        if var in temporal_sources:
            check2 = [startdate >= valid_dates[source][0] for source in sources]
            check3 = [enddate <= valid_dates[source][1] for source in sources]
        else:
            check2 = [True]
            check3 = [True]

        results[var] = {source: {"valid_source:": check1[i],
                                 "valid_startdate:": check2[i],
                                 "valid_enddate:": check3[i],
                                } for i, source in enumerate(sources)}

        all_results.append(np.all([check1, check2, check3]))
    
    succes = np.all(all_results)

    return results, succes

def nasa_account(user_pw = None):

    n_max = 3
    succes = False
    n = 1

    while not succes and n <= n_max:

        folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
        test_url = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NXGLC.5.12.4/1987/08/MERRA2_100.tavg3_2d_glc_Nx.19870801.nc4"
        test_file = os.path.join(folder, "nasa_test.nc4")

        if os.path.isfile(test_file):
            os.remove(test_file)

        if not isinstance(user_pw, type(None)):
            username, password = user_pw
        else:
            username, password = pywapor.collect.get_pw_un.get('NASA')

        x = requests.get(test_url, allow_redirects = False)
        y = requests.get(x.headers['location'], auth = (username, password))

        if x.ok and y.ok:

            with open(test_file, 'w+b') as z:
                z.write(y.content)
        
            if os.path.isfile(test_file):
                statinfo = os.stat(test_file)
                succes = statinfo.st_size == 3963517
                os.remove(test_file)
                if not succes:
                    error = "please add 'NASA GESDISC DATA ARCHIVE' to 'Approved Applications'."
            else:
                succes = False

        else:
            error = "wrong username/password."
            succes = False
            if os.path.isfile(test_file):
                os.remove(test_file)

        if not succes:
            print(f"Try {n}/{n_max} failed: {error}")
            time.sleep(3)
            
        n += 1

    return succes

# nasa_account()