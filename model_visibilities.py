# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 09:11:22 2021

@author: jdrevon
"""

import numpy as np
from A1_mas_to_rad import mas_to_rad
from scipy.special import j1,jv

def V_uniform(q,diam):
    
    V = 2*j1(np.pi*q*mas_to_rad(diam))/(np.pi*q*mas_to_rad(diam))
    
    return V

def V_ring(q,inner_diam, outter_diam):
    
    if inner_diam == 0:    
        V = V_uniform(q,outter_diam)
    else:            
        f = ((outter_diam-inner_diam))/(inner_diam)
        factor = 2/(np.pi*q*mas_to_rad(inner_diam)*(2*f+f**2))
        V = factor*((1+f)*j1((1+f)*mas_to_rad(inner_diam)*np.pi*q)-j1(mas_to_rad(inner_diam)*np.pi*q))
        
    return V


def V2_uniform(q,diam):
    
    V2 = V_uniform(q,diam)**2
    
    return V2

def V_gaussian(q,FWHM):

    V = np.exp((-(np.pi*mas_to_rad(FWHM)*q)**2)/(4*np.log(2)))
    
    return V


def V2_gaussian(q,FWHM):

    V2 = V_gaussian(q,FWHM)**2
    
    return V2

def V_LD(q,diam,u_lambda):
    
    x = np.pi*q*mas_to_rad(diam)
    alpha = 1-u_lambda 
    beta = u_lambda
    V = (alpha*j1(x)/x+beta*(np.pi/2)**(1/2)*jv(3/2,x)/x**(3/2))/(alpha/2+beta/3)
    
    return V


# plt.figure()
# plt.scatter(q_UD,V_LD(q_UD,10,0.5)**2,s=1)
# plt.scatter(q_UD,V_uniform(q_UD,10)**2,s=1)
# ax=plt.gca()
# ax.set_xscale('log')
# ax.set_yscale('log')