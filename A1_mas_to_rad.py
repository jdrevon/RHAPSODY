#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 12:41:29 2020

@author: jdrevon
"""
import numpy as np

def mas_to_rad(diam_star):
    y = diam_star*1E-3*np.pi/(180*3600)
    return y 

def rad_to_mas(rad):
    y = rad*1E3/(np.pi/(180*3600))
    return y 

def au_to_R_sun(star_radius):
    
    R_final =  star_radius*214.9394693836
    
    return R_final