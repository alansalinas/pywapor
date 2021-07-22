import requests
import watertools
import os
import time

def nasa_account(user_pw = None):

    n_max = 3
    succes = False
    n = 1

    while not succes and n <= n_max:

        test_url = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NXGLC.5.12.4/1987/08/MERRA2_100.tavg3_2d_glc_Nx.19870801.nc4"
        test_file = os.path.join(os.getcwd(), "nasa_test.nc4")

        if os.path.isfile(test_file):
            os.remove(test_file)

        if not isinstance(user_pw, type(None)):
            username, password = user_pw
        else:
            username, password = watertools.Functions.Random.Get_Username_PWD.GET('NASA')

        x = requests.get(test_url, allow_redirects = False)
        y = requests.get(x.headers['location'], auth = (username, password))

        if x.ok and y.ok:

            with open(test_file, 'wb') as z:
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

        if not succes:
            print(f"Try {n}/{n_max} failed: {error}")
            time.sleep(3)
            
        n += 1

    return succes