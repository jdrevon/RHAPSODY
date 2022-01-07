#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 19:19:11 2021

@author: jdrevon
"""

import numpy as np
from brightness_distribution_models import normalized_brightness_R

def check_intensity(res_F_mod, res_UD, res_UD_bef, intensity_init, flux_init):

        I_tot_tmp = (res_F_mod*np.array([normalized_brightness_R(res_UD[i]/2, res_UD_bef[i], res_UD[i])[0] for i in range(len(res_UD))]))
        I_tot_norm = I_tot_tmp

        r_prob_index = np.where(I_tot_norm>1)[0]
                  
        F_new = res_F_mod
                
        if np.shape(r_prob_index)[0]!= 0:
            for i in range(0,len(r_prob_index)):
                index = r_prob_index[i]
                F_new[index] = flux_init[index]
                I_tot_norm[index] =  intensity_init[index]
            
        return F_new, I_tot_norm