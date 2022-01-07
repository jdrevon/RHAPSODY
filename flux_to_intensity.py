# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:45:39 2022

@author: jdrevon
"""

import numpy as np
from brightness_distribution_models import normalized_brightness_R

def flux_to_intensity(min_ring,max_ring,flux):
    

    intensity = (flux*np.array([normalized_brightness_R(max_ring[i]/2, min_ring[i], max_ring[i])[0] for i in range(len(min_ring))]))
    
    return intensity
    


    