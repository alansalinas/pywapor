# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:31:37 2020

@author: timhe
"""
import os
import json
import pywapor
from cryptography.fernet import Fernet

def get(server):

    folder = os.path.dirname(os.path.realpath(pywapor.__path__[0]))
    filename = "secret.txt"
    key_file = os.path.join(folder, filename)
    json_file = os.path.join(folder, "keys.json")
    
    if not os.path.exists(key_file):
        pywapor.collect.setup_dl_accounts.create_key()
        if os.path.exists(json_file):
            print("removing old/invalid json file")
            os.remove(json_file)
    
    if not os.path.exists(json_file):
        pywapor.collect.setup_dl_accounts.setup_account(server)

    f = open(key_file,"r")
    key = f.read()
    f.close()
    
    cipher_suite = Fernet(key.encode('utf-8'))    
    
    # open json file  
    with open(json_file) as f:
        datastore = f.read()
    obj = json.loads(datastore)      
    f.close()

    if server not in obj.keys():
        pywapor.collect.setup_dl_accounts.setup_account(server)
        obj = None 
        with open(json_file) as f:
            datastore = f.read()
        obj = json.loads(datastore)      
        f.close()       
    
    if not server == "WAPOR":
        username_crypt, pwd_crypt = obj[server]
        username = cipher_suite.decrypt(username_crypt.encode("utf-8"))
        pwd = cipher_suite.decrypt(pwd_crypt.encode("utf-8"))        
    else:
        username_crypt = obj[server]
        username = cipher_suite.decrypt(username_crypt[0].encode("utf-8"))
        pwd = b''
    
    return(str(username.decode("utf-8")), str(pwd.decode("utf-8")))