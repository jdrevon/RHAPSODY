# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:51:38 2021

@author: jdrevon
"""

import numpy as np
from pretty_table_reading import READ_PRETTY_TABLE

def output_intensity_profile(PATH_INTENSITY_LM, PATH_INTENSITY_N, wavel_UD):

    header_LM, data_LM = READ_PRETTY_TABLE(PATH_INTENSITY_LM,1)
    header_N, data_N   = READ_PRETTY_TABLE(PATH_INTENSITY_N,1)

    wavel_LM =  header_LM[1:].astype('float')
    wavel_N  =  header_N[1:].astype('float')

    distance_LM = data_LM[:,0]
    distance_N  = data_N[:,0]
    
    intensity_LM = data_LM[:,1:]    
    intensity_N  = data_N[:,1:] 
    
    return distance_LM, intensity_LM, distance_N, intensity_N, wavel_LM, wavel_N