#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:19:23 2021

@author: jdrevon
"""

import numpy as np
from brightness_distribution_models import normalized_brightness_R
from A2_gaussian import gaussian
from rhapsody_init import sigma, init, alpha, DATA_band_name, size_ring, NBR_ring


import sys 

def rings():

    """

    Outputs: O
    
    diam_UD_bef [mas] : Inner diameter of the generated rings for each bands
    diam_UD [mas] : Outter diameter of the generated rings for each bands
    
    I_ratio : Guess profile of the radial intensity in each bands
      
    flux_ratio [mas] : Flux associated to each rings based on the guess profile of the radial intensity in each band        
     
    """    
    
    # 1 Initialization of the ring size with respect to the user's parameters
    
    #I create an array to set the outter boundaries of each rings
    
    diam_UD = [np.round(np.linspace(size_ring[i],size_ring[i]*NBR_ring[i], NBR_ring[i]),2) for i in range(len(DATA_band_name))]
    
    # I create an array to set the inner boundaries of each rings

    diam_UD_bef= [np.array([0]+diam_UD[i].tolist())[:-1] for i in range(len(diam_UD))]
    
    # 2 Initialization of the initial intensity profile with respect to the power law given by the user with alpha    
  
    if init == 'PW':  
  
        I_ratio = [1/diam_UD[i]**alpha[i]/(1/diam_UD[i][0]**alpha[i]) for i in range(len(diam_UD))]

    elif init == 'G':  
        
        I_ratio = [gaussian(0, sigma[i], 1, diam_UD[i])/gaussian(0,sigma[i], 1, diam_UD[i][0]) for i in range(len(diam_UD))]
        
    else: 
        print("Please enter a valid value: Power (PW) or Gaussian (G) Law")
        sys.exit()

    # 3 Determination of the associated flux of each rings with respect to the user's intensity profile and rings' dimensions

    flux = [I_ratio[k]/np.array([normalized_brightness_R(diam_UD[k][i]/2, diam_UD_bef[k][i], diam_UD[k][i])[0] for i in range(len(diam_UD[k]))]) for k in range(len(diam_UD))]
    flux_max = [1/np.array([normalized_brightness_R(diam_UD[k][i]/2, diam_UD_bef[k][i], diam_UD[k][i])[0] for i in range(len(diam_UD[k]))]) for k in range(len(diam_UD))]         

    return diam_UD_bef, diam_UD, I_ratio, flux, flux_max
