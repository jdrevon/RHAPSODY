#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 09:39:21 2021

@author: jdrevon
"""
import numpy as np
from scipy.constants import c,h,k
from astropy.constants import R_sun

def BB_v(wavel,T):
    
    freq = np.array(c/wavel).astype('float64')
    black_body = 2*h*freq**3/c**2*1/(np.exp(h*freq/(k*T))-1)
    return black_body #*freq*(np.pi*mas_to_rad(10.15)**2/4) #factor for flux


def BB_v_Jy(wavel, T_eff, R_sun_usr, distance):
    
    black_body = BB_v(wavel,T_eff)*4*np.pi*(R_sun_usr*R_sun.value)**2*np.pi/(1E-26)/(4*np.pi*(distance*3.08E16)**2)
    
    return wavel, black_body



# wavel = np.linspace(0.1,20, 1000)
# T_eff = 1E4


# # BB = BB_v_Jy(wavel, T_eff, R_sun_usr, distance)
# import matplotlib.pyplot as plt

# BB = BB_v(wavel*1E-6, T_eff)
# plt.figure()
# plt.plot(wavel, BB)
# plt.yscale('log')
# plt.xscale('log')