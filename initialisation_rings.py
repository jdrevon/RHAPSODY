#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 15:19:23 2021

@author: jdrevon
"""

import numpy as np
from brightness_distribution_models import normalized_brightness_R
from A2_gaussian import gaussian
from rhapsody_init import sigma_LM, sigma_N, init

import sys 

def rings(size_ring_LM, size_ring_N, NBR_ring_LM, NBR_ring_N, NBR_WL_LM, NBR_WL_N, alpha_LM, alpha_N):

    """
    
    Inputs: I
    
    size_ring_LM [mas] : constant width of a ring (outter diameter-inner diameter) in LM band
    size_ring_LM [mas] : constant width of a ring (outter diameter-inner diameter) in N band
    
    NBR_ring_LM [#] : number of rings to put in the model in LM band
    NBR_ring_N  [#] : number of rings to put in the model in N band
    
    NBR_WL_LM [#] : number of wavelengths to model in LM band
    NBR_WL_N  [#] : number of wavelengths to model in N band
                
    alpha_LM [#] : Power law of the intensity profile for LM band (1/r**(alpha))
    alpha_N  [#] : Power law of the initial guess intensity profile for N band (1/r**(alpha))

    Outputs: O
    
    diam_UD_bef [mas] : Inner diameter of the generated rings for LM band (diam_UD[0]) and N band (diam_UD[1])
    diam_UD [mas] : Outter diameter of the generated rings for LM band (diam_UD[0]) and N band (diam_UD[1])
    
    I_ratio_LM [mas] : Guess profile of the radial intensity in LM band
    I_ratio_N [mas] : Guess profile of the radial intensity in N band
      
    flux_ratio_LM [mas] : Flux associated to each rings based on the guess profile of the radial intensity in LM band        
    flux_ratio_N  [mas] : Flux associated to each rings based on the guess profile of the radial intensity in N band   
     
    """    
    
    # 1 Initialization of the ring size with respect to the user's parameters
      
    
    #I create an array to set the outter boundaries of each rings
    
    diam_UD_LM = np.round(np.linspace(size_ring_LM,size_ring_LM*NBR_ring_LM, NBR_ring_LM),2)
    
    diam_UD_N = np.round(np.linspace(size_ring_N,size_ring_N*NBR_ring_N,NBR_ring_N),2)

    diam_UD = [diam_UD_LM,diam_UD_N]

    # I create an array to set the inner boundaries of each rings

    diam_UD_bef= np.array([np.array([0]+diam_UD_LM.tolist())[:-1],np.array([0]+diam_UD_N.tolist())[:-1]],dtype='object')
    
    diam_UD_bef_LM = diam_UD_bef[0]
    diam_UD_bef_N = diam_UD_bef[1]

    # 2 Initialization of the initial intensity profile with respect to the power law given by the user with alpha    
  
    if init == 'PW':  
  
        I_ratio_LM = 1/diam_UD_LM**alpha_LM/(1/diam_UD_LM[0]**alpha_LM)
        I_ratio_N  = 1/diam_UD_N**alpha_N/(1/diam_UD_N[0]**alpha_N)

    elif init == 'G':  
        
        I_ratio_LM = gaussian(0, sigma_LM, 1, diam_UD_LM)/gaussian(0,sigma_LM, 1, diam_UD_LM[0])
        I_ratio_N  = gaussian(0, sigma_N, 1, diam_UD_N)/gaussian(0, sigma_N, 1, diam_UD_N[0]) 
        
    else: 
        print("Please enter a valid value: Power (PW) or Gaussian (G) Law")
        sys.exit()

    # 3 Determination of the associated flux of each rings with respect to the user's intensity profile and rings' dimensions

    flux_LM  = I_ratio_LM/np.array([normalized_brightness_R(diam_UD_LM[i]/2, diam_UD_bef_LM[i], diam_UD_LM[i])[0] for i in range(len(diam_UD_LM))])
    flux_N   = I_ratio_N/np.array([normalized_brightness_R(diam_UD_N[i]/2, diam_UD_bef_N[i], diam_UD_N[i])[0] for i in range(len(diam_UD_N))])
             

    flux_LM_max = 1/np.array([normalized_brightness_R(diam_UD_LM[i]/2, diam_UD_bef_LM[i], diam_UD_LM[i])[0] for i in range(len(diam_UD_LM))])
    flux_N_max  = 1/np.array([normalized_brightness_R(diam_UD_N[i]/2, diam_UD_bef_N[i], diam_UD_N[i])[0] for i in range(len(diam_UD_N))])

    # flux_LM = flux_LM/np.sum(flux_LM)
    # flux_N  = flux_N/np.sum(flux_N)
    
    # flux_LM_max = flux_LM_max/np.sum(flux_LM)
    # flux_N_max = flux_N_max/np.sum(flux_N)

    return diam_UD_bef, diam_UD, I_ratio_LM, I_ratio_N, flux_LM, flux_N, flux_LM_max, flux_N_max
