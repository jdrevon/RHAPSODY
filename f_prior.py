# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 19:39:49 2021

@author: jdrevon
"""

import numpy as np

def total_variation(hp,intensity,radius):

    delta_I  = np.diff(intensity) 
    delta_r  = np.diff(radius)
    gradient_2 = np.abs(delta_I/delta_r)**2
    f_prior = hp*np.sum(gradient_2)
    
    return f_prior
    
    
def quad_smoothness(hp,intensity):
    
    return hp*np.sum(np.sqrt(np.diff(intensity)**2))