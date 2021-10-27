#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 11:30:25 2020

@author: jdrevon
"""
import numpy as np
import numpy.ma as ma
 
def nearest_index(array, value):
    
    array = np.array(array)
    idx = np.abs(array - value).argmin()
    return idx

    

# before,after=nearest_index(wavel,lambda_tot[0])
# print(lambda_tot[0])
# print(wavel[before][0])
# print(wavel[after][0])
