import os
import json
import pywapor
import getpass
import sys
import requests
import time
from pywapor.general.logger import log
from cryptography.fernet import Fernet

def setup(account):
    """Asks, saves and tests a username/password combination for `account`.

    Parameters
    ----------
    account : {"NASA" | "VITO" | "WAPOR"}
        Which un/pw combination to store.
    """

    folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
    filename = "secret.txt"
    filehandle = os.path.join(folder, filename)
    json_filehandle = os.path.join(folder, "keys.json")

    if not os.path.exists(filehandle):
        create_key()
    
    f = open(filehandle,"r")
    key = f.read()
    f.close()
    
    cipher_suite = Fernet(key.encode('utf-8'))

    if os.path.exists(json_filehandle):
        # open json file  
        with open(json_filehandle) as f:
            datastore = f.read()
        obj = json.loads(datastore)      
        f.close()
    else:
        obj = {}

    if account == "WAPOR":
        API_Key = input(f"{account} API token: ")
        API_Key_crypt = cipher_suite.encrypt(("%s" %API_Key).encode('utf-8'))
        obj[account] = [str(API_Key_crypt.decode("utf-8"))]
    else:
        account_name = input(f"{account} username: ")
        pwd = getpass.getpass(f"{account} password: ")            
        account_name_crypt = cipher_suite.encrypt(("%s" %account_name).encode('utf-8'))
        pwd_crypt = cipher_suite.encrypt(("%s" %pwd).encode('utf-8'))
        obj[account] = ([str(account_name_crypt.decode("utf-8")), str(pwd_crypt.decode("utf-8"))])

    # save extent in task
    with open(json_filehandle, 'w') as outfile:
        json.dump(obj, outfile)        

    if account == "NASA":
        log.info("--> Testing NASA un/pw.")
        succes = nasa_account()
        if succes:
            log.info("--> NASA un/pw working.")
        else:
            _ = obj.pop("NASA")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your NASA account.")            

    if account == "VITO":
        log.info("--> Testing VITO un/pw.")
        succes = vito_account()
        if succes:
            log.info("--> VITO un/pw working.")
        else:
            _ = obj.pop("VITO")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your VITO account.")  

    if account == "WAPOR":
        log.info("--> Testing WAPOR token.")
        succes = wapor_account()
        if succes:
            log.info("--> WAPOR token working.")
        else:
            _ = obj.pop("WAPOR")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your WAPOR token.")  

    return

def get(account):
    """Loads a required username/password.

    Parameters
    ----------
    account : {"NASA" | "VITO" | "WAPOR"}
        Which un/pw combination to load.
    """

    folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
    filename = "secret.txt"
    key_file = os.path.join(folder, filename)
    json_file = os.path.join(folder, "keys.json")
    
    if not os.path.exists(key_file):
        create_key()
        if os.path.exists(json_file):
            print("removing old/invalid json file")
            os.remove(json_file)
    
    if not os.path.exists(json_file):
        setup(account)

    f = open(key_file,"r")
    key = f.read()
    f.close()
    
    cipher_suite = Fernet(key.encode('utf-8'))    
    
    # open json file  
    with open(json_file) as f:
        datastore = f.read()
    obj = json.loads(datastore)      
    f.close()

    if account not in obj.keys():
        setup(account)
        obj = None 
        with open(json_file) as f:
            datastore = f.read()
        obj = json.loads(datastore)      
        f.close()       
    
    if not account == "WAPOR":
        username_crypt, pwd_crypt = obj[account]
        username = cipher_suite.decrypt(username_crypt.encode("utf-8"))
        pwd = cipher_suite.decrypt(pwd_crypt.encode("utf-8"))        
    else:
        username_crypt = obj[account]
        username = cipher_suite.decrypt(username_crypt[0].encode("utf-8"))
        pwd = b''
    
    return(str(username.decode("utf-8")), str(pwd.decode("utf-8")))

def create_key():
    """Generates a key file.
    """
    folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
    filename = "secret.txt"
    filehandle = os.path.join(folder, filename)

    if os.path.exists(filehandle):
        os.remove(filehandle)
    
    f = open(filehandle,"w+")
    f.write(str(Fernet.generate_key().decode("utf-8")))
    f.close()

def vito_account(user_pw = None):
    """Check if the given or stored VITO username/password combination
    is correct. Accounts can be created on https://www.vito-eodata.be.

    Parameters
    ----------
    user_pw : tuple, optional
        ("username", "password") to check, if `None` will try to load the 
        password from the keychain, by default None.

    Returns
    -------
    bool
        True if the password works, otherwise False.
    """
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
            username, password = get('VITO')

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
    """Check if the given or stored NASA username/password combination is 
    correct. Accounts can be created on https://urs.earthdata.nasa.gov/users/new.

    Parameters
    ----------
    user_pw : tuple, optional
        ("username", "password") to check, if `None` will try to load the 
        password from the keychain, by default None.

    Returns
    -------
    bool
        True if the password works, otherwise False.
    """
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
            username, password = get('NASA')

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
    """Check if the given or stored WAPOR token is 
    correct. Accounts can be created on https://wapor.apps.fao.org/home/WAPOR_2/1.

    Parameters
    ----------
    user_pw : tuple, optional
        ("", "token") to check, if `None` will try to load the 
        password from the keychain, by default None.

    Returns
    -------
    bool
        True if the password works, otherwise False.
    """
    n_max = 3
    succes = False
    n = 1

    while not succes and n <= n_max:

        if not isinstance(user_pw, type(None)):
            username, _ = user_pw
        else:
            username, _ = get('WAPOR')

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