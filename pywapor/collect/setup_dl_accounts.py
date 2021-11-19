# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 08:28:08 2020

@author: timhe
"""
import os
import json
from cryptography.fernet import Fernet
import pywapor
import getpass
import sys

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
        succes = pywapor.general.tests.nasa_account()
        if succes:
            print("NASA account: working")
        else:
            to_be_removed = obj.pop("NASA")
            with open(json_filehandle, 'w') as outfile:
                json.dump(obj, outfile)
            sys.exit(f"Please fix your NASA account.")            
   
    return
    
