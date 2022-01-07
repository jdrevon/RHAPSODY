#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 22:29:05 2021

@author: jdrevon
"""

import os
import glob 
import shutil

def COPY_DATA(BEFORE, AFTER):
    
    if not os.path.exists(AFTER):
        os.makedirs(AFTER)

    if os.path.exists(AFTER+'/LM'):
        for filenames in glob.glob(AFTER+'/LM/*.fits'):
            os.remove(filenames)
            

    if os.path.exists(AFTER+'/N'):
        for filenames in glob.glob(AFTER+'/N/*.fits'):
            os.remove(filenames)
    
    if os.path.exists(AFTER+'/TRASH'):
        for filenames in glob.glob(AFTER+'/TRASH/*.fits'):
            os.remove(filenames)


    for filenames in glob.glob(BEFORE+'/*.fits'):
        shutil.copy2(filenames,AFTER+'/')

    return 