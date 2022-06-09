#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:03:17 2021

@author: jdrevon
"""
import numpy as np
from pretty_table_reading import READ_PRETTY_TABLE
from brightness_distribution_models import normalized_brightness_R, normalized_brightness_UD
from scipy.interpolate import interp2d, interp1d
def increase_intensity_profile_res(PATH_INTENSITY_PROFILE,**kwargs):

    
        header, data = READ_PRETTY_TABLE(PATH_INTENSITY_PROFILE,1)
        
        wavel = (header[1:]).astype(np.float)
        inner_radius = (data[:,0]-np.diff(data[:,0])[0])
        intensity = data[:,1:].T
        
        min_dist        = kwargs.get('min_dist', min(inner_radius))    
        max_dist        = kwargs.get('max_dist', max(inner_radius))    
        nb_points       = kwargs.get('nb_points', len(inner_radius))    

        rho   = np.linspace(min_dist,max_dist,nb_points)

        function_I = interp2d(inner_radius,wavel,intensity)
        
        I_tot_norm_HR  = function_I(rho, wavel)

                            
        return rho, I_tot_norm_HR, wavel


def increase_intensity_profile_res_1W(PATH_INTENSITY_PROFILE,**kwargs):

    
        header, data = READ_PRETTY_TABLE(PATH_INTENSITY_PROFILE,1)
        
        wavel = (header[1:]).astype(np.float)
        inner_radius = (data[:,0]-np.diff(data[:,0])[0])
        intensity = data[:,1:].T
        
        min_dist        = kwargs.get('min_dist', min(inner_radius))    
        max_dist        = kwargs.get('max_dist', max(inner_radius))    
        nb_points       = kwargs.get('nb_points', len(inner_radius))    

        rho   = np.linspace(min_dist,max_dist,nb_points)

        function_I = interp1d(inner_radius,intensity)
        
        I_tot_norm_HR  = function_I(rho)

                            
        return rho, I_tot_norm_HR, wavel
