import requests
import pywapor
import os
import time
from pywapor.general.logger import log

def password_reqs():
    password_reqs = {"NASA": ["MOD13", "MYD13", "MCD43", 
                                "MOD11", "MYD11", "MERRA2"],
                    "VITO": ["PROBAV"],
                    "WAPOR": ["WAPOR"],}
    return password_reqs

def vito_account(user_pw = None):

    n_max = 3
    succes = False
    n = 1

    while not succes and n <= n_max:

        folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
        test_url = r"https://www.vito-eodata.be/PDF/datapool/Free_Data/PROBA-V_300m/S1_TOC_-_300_m_C1/2014/1/15/PV_S1_TOC-20140115_333M_V101/PROBAV_S1_TOC_20140115_333M_V101.VRG"
        test_file = os.path.join(folder, "vito_test.vrg")

        if os.path.isfile(test_file):
            os.remove(test_file)

        if not isinstance(user_pw, type(None)):
            username, password = user_pw
        else:
            username, password = pywapor.collect.get_pw_un.get('VITO')

        x = requests.get(test_url, auth = (username, password))

        if x.ok:

            with open(test_file, 'w+b') as z:
                z.write(x.content)

                statinfo = os.stat(test_file)
                succes = statinfo.st_size == 15392
                os.remove(test_file)
                if not succes:
                    error = "something went wrong."

        else:
            error = "wrong username/password."
            succes = False
            if os.path.isfile(test_file):
                os.remove(test_file)

        if not succes:
            log.warning(f"Try {n}/{n_max} failed.")
            time.sleep(3)
            
        n += 1

    if not succes:
        log.warning(error)

    return succes

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
            log.warning(f"Try {n}/{n_max} failed.")
            time.sleep(3)
            
        n += 1

    if not succes:
        log.warning(error)

    return succes

def wapor_account(user_pw = None):

    n_max = 3
    succes = False
    n = 1

    while not succes and n <= n_max:

        if not isinstance(user_pw, type(None)):
            username, _ = user_pw
        else:
            username, _ = pywapor.collect.get_pw_un.get('WAPOR')

        sign_in= 'https://io.apps.fao.org/gismgr/api/v1/iam/sign-in'
        resp_vp=requests.post(sign_in,headers={'X-GISMGR-API-KEY':username})
        
        if resp_vp.ok:
            succes = True
        else:
            error = "wrong token."
            succes = False

        if not succes:
            log.warning(f"Try {n}/{n_max} failed.")
            time.sleep(3)
            
        n += 1

    if not succes:
        log.warning(error)

    return succes

if __name__ == "__main__":
    nasa_succes = nasa_account()
    vito_succes = vito_account()
    wapor_succes = wapor_account()

    print(nasa_succes, vito_succes, wapor_succes)