import os
import json
from cryptography.fernet import Fernet
import pywapor
import getpass
import sys
from pywapor.general.logger import log

def create_key():

    folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
    filename = "secret.txt"
    filehandle = os.path.join(folder, filename)

    if os.path.exists(filehandle):
        os.remove(filehandle)
    
    f = open(filehandle,"w+")
    f.write(str(Fernet.generate_key().decode("utf-8")))
    f.close()

    return

def setup_account(account):

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
        succes = pywapor.general.tests.nasa_account()
        if succes:
            log.info("--> NASA un/pw working.")
        else:
            _ = obj.pop("NASA")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your NASA account.")            

    if account == "VITO":
        log.info("--> Testing VITO un/pw.")
        succes = pywapor.general.tests.vito_account()
        if succes:
            log.info("--> VITO un/pw working.")
        else:
            _ = obj.pop("VITO")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your VITO account.")  

    if account == "WAPOR":
        log.info("--> Testing WAPOR token.")
        succes = pywapor.general.tests.wapor_account()
        if succes:
            log.info("--> WAPOR token working.")
        else:
            _ = obj.pop("WAPOR")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your WAPOR token.")  

    return