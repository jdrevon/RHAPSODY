#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 15:04:53 2021

@author: jdrevon
"""


import glob,os
import shutil

def REMOVE_ALL_FILES_FROM_FOLDER(FOLDER):

    try :
        files = glob.glob(FOLDER+'*')
    except : 
        files = glob.glob(FOLDER+'/*')
        
    for f in files:
        try :
            os.remove(f)
        except : 
            shutil.rmtree(f, ignore_errors=True)
    return 
