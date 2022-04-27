#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 22:29:05 2021

@author: jdrevon
"""

import os
import glob 
import shutil
from rhapsody_init import DATA_band_name

def COPY_DATA(BEFORE, AFTER):
    
    print('CREATING NEW DIRECTORIES')
    
    if not os.path.exists(AFTER):
        os.makedirs(AFTER)

    for i in range(len(DATA_band_name)):
                
        if os.path.exists(AFTER+'/'+DATA_band_name[i]):
            for filenames in glob.glob(AFTER+'/'+DATA_band_name[i]+'/*.fits'):
                os.remove(filenames)
        else: 
            os.makedirs(AFTER+'/'+DATA_band_name[i])

                
        if os.path.exists(AFTER+'/TRASH'):
            for filenames in glob.glob(AFTER+'/TRASH/*.fits'):
                os.remove(filenames)
       
        else:
            os.makedirs(AFTER+'/TRASH')
    

        for filenames in glob.glob(BEFORE[i]+'/*.fits'):
            shutil.copy2(filenames,AFTER+'/'+DATA_band_name[i]+'/')

    return 