# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 17:53:21 2021

@author: jdrevon
"""

import os

def FOLDER_CREATION(PATH_OUTPUT):
    
    if not os.path.exists(PATH_OUTPUT):
        os.makedirs(PATH_OUTPUT)

    return 