# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 16:38:30 2021

@author: jdrevon
"""

import numpy as np
from A6_nearest_index import nearest_index

def normalized_brightness_UD(rho, diam):
    if np.logical_or(isinstance(rho,float),isinstance(rho,int)):
        rho = np.array([rho])
    I_UD = np.zeros(len(rho)) 
    I_UD[rho<=diam/2] += 4/(np.pi*diam**2)
    
    return I_UD

def normalized_brightness_R(rho, inner_diam, outter_diam):
    if np.logical_or(isinstance(rho,float),isinstance(rho,int)):
        rho = np.array([rho])
    I_UD = np.zeros(len(rho)) 
    
    if inner_diam == 0:
        I_UD = normalized_brightness_UD(rho,outter_diam)        
    
    else : 
        cond  = ((rho>inner_diam/2) & (rho<=outter_diam/2)) # fonction porte
        I_UD[cond] += 4/(np.pi*((outter_diam)**2-(inner_diam)**2))
    
    return I_UD


def normalized_brightness_Gaussian(rho, FWHM):
    
    I_G = 4*np.log(2)/(np.pi*FWHM**2)*np.exp(-(rho/FWHM)**2*4*np.log(2))
    
    return I_G

def normalized_brightness_LD(rho, diam, u_lambda):
    
    cond = rho<diam/2 # fonction porte
    
    I_LD = np.array([1E-10]*(len(rho))) 
    
    alpha  =1-4/np.pi+8/(np.pi**2)
    I_LD[cond] =  4/(np.pi*diam**2)*1/(1-alpha*u_lambda)*(1-u_lambda*(1-np.cos(np.pi*rho[cond]/diam)))
    
    return I_LD


def intensity_1UD(rho, diam):
    
    I_UD = normalized_brightness_UD(rho,diam)
    
    return rho, I_UD

def intensity_2UD(rho, diam, diam2, flux_ratio_UD2, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]
    I_UD  = normalized_brightness_UD(rho,diam)
    I_UD2 = normalized_brightness_UD(rho,diam2)
    
    I_tot = F_tot/(1+flux_ratio_UD2)*I_UD[0]*(I_UD/I_UD[0]+flux_ratio_UD2*I_UD2/I_UD[0])  
    
    return rho,I_tot, I_UD,I_UD2


def intensity_3UD(rho, diam, diam2, diam3, flux_ratio_UD2, flux_ratio_UD3, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]
    I_UD  = normalized_brightness_UD(rho,diam)
    I_UD2 = normalized_brightness_UD(rho,diam2)
    I_UD3 = normalized_brightness_UD(rho,diam3)
    
    I_tot = F_tot/(1+flux_ratio_UD2+ flux_ratio_UD3)*I_UD[0]*(I_UD/I_UD[0]+flux_ratio_UD2*I_UD2/I_UD[0]+flux_ratio_UD3*I_UD3/I_UD[0])  
    
    return rho, I_tot

def intensity_4UD(rho, diam, diam2, diam3, diam4, flux_ratio_UD2, flux_ratio_UD3, flux_ratio_UD4, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]
    I_UD  = normalized_brightness_UD(rho,diam)
    I_UD2 = normalized_brightness_UD(rho,diam2)
    I_UD3 = normalized_brightness_UD(rho,diam3)
    I_UD4 = normalized_brightness_UD(rho,diam4)
    
    I_tot = F_tot/(1+flux_ratio_UD2+ flux_ratio_UD3 + flux_ratio_UD4)*I_UD*(1+flux_ratio_UD2*I_UD2/I_UD+flux_ratio_UD3*I_UD3/I_UD+flux_ratio_UD4*I_UD4/I_UD)  
    
    return rho, I_tot

def intensity_6UD(rho, diam, diam2, diam3, diam4, diam5, diam6, flux_ratio_UD2, flux_ratio_UD3, flux_ratio_UD4, flux_ratio_UD5, flux_ratio_UD6, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]
    I_UD  = normalized_brightness_UD(rho,diam)
    I_UD2 = normalized_brightness_UD(rho,diam2)
    I_UD3 = normalized_brightness_UD(rho,diam3)
    I_UD4 = normalized_brightness_UD(rho,diam4)
    I_UD5 = normalized_brightness_UD(rho,diam5)
    I_UD6 = normalized_brightness_UD(rho,diam6)
    
    I_tot = F_tot/(1+flux_ratio_UD2+ flux_ratio_UD3 + flux_ratio_UD4 + flux_ratio_UD5 + flux_ratio_UD6)\
        *I_UD*(1+flux_ratio_UD2*I_UD2/I_UD+flux_ratio_UD3*I_UD3/I_UD+flux_ratio_UD4*I_UD4/I_UD+flux_ratio_UD5*I_UD5/I_UD+flux_ratio_UD6*I_UD6/I_UD)  
    
    return rho, I_tot

def intensity_1GAUSSIAN(rho, FWHM):
    
    I_tot = normalized_brightness_Gaussian(rho,FWHM)
    
    return rho, I_tot

def intensity_1UD_1GAUSSIAN(rho, diam, FWHM, flux_ratio_G1, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]

    I_UD  = normalized_brightness_UD(rho,diam)
    I_G1   = normalized_brightness_Gaussian(rho,FWHM)
    
    I_tot = F_tot/(1+flux_ratio_G1)*I_UD*(1+flux_ratio_G1*I_G1/I_UD)  
    
    return rho, I_tot


def intensity_2UD_1GAUSSIAN(rho, diam, diam2, FWHM, flux_ratio_UD2, flux_ratio_G1, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]

    I_UD   = normalized_brightness_UD(rho,diam)
    I_UD2  = normalized_brightness_UD(rho,diam2)
    I_G1   = normalized_brightness_Gaussian(rho,FWHM)
    
    I_tot = F_tot/(1+flux_ratio_G1+ flux_ratio_UD2)*I_UD*(1+flux_ratio_G1*I_G1/I_UD+flux_ratio_UD2*I_UD2/I_UD)  
    
    return rho, I_tot


def intensity_1LD(rho, diam_LD, u_lambda):
    
    I_tot   = normalized_brightness_LD(rho, diam_LD, u_lambda)
    
    return rho, I_tot

def intensity_1LD_1UD(rho, diam_UD, diam_LD, flux_ratio_LD, u_lambda, ISO_spectra_x, ISO_spectra_y, lambda_c):

    index = nearest_index(ISO_spectra_x,lambda_c)
    F_tot = ISO_spectra_y[index]
    
    I_LD   = normalized_brightness_LD(rho, diam_LD, u_lambda)
    I_UD   = normalized_brightness_UD(rho, diam_UD)
    
    I_tot = F_tot/(1+flux_ratio_LD)*I_UD*(1+flux_ratio_LD*I_LD/I_UD)  
    
    return rho, I_tot

def intensity_2UD_1LD(rho, diam_UD1, diam_UD2, diam_LD, u_lambda, flux_ratio_UD2, flux_ratio_LD, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x, lambda_c)
    F_tot = ISO_spectra_y[index]

    I_UD   = normalized_brightness_UD(rho, diam_UD1)
    I_UD2  = normalized_brightness_UD(rho, diam_UD2)
    I_LD   = normalized_brightness_LD(rho, diam_LD, u_lambda)
    
    I_tot = F_tot/(1+flux_ratio_UD2+ flux_ratio_LD)*I_UD*(1+flux_ratio_UD2*I_UD2/I_UD+flux_ratio_LD*I_LD/I_UD)  
    
    return rho, I_tot

def intensity_1UD_1LD_1GAUSSIAN(rho, diam_UD1, diam_G1, diam_LD1, u_lambda, flux_ratio_G1, flux_ratio_LD1, ISO_spectra_x, ISO_spectra_y, lambda_c):
    
    index = nearest_index(ISO_spectra_x, lambda_c)
    F_tot = ISO_spectra_y[index]

    I_UD   = normalized_brightness_UD(rho, diam_UD1)
    I_G1   = normalized_brightness_UD(rho, diam_G1)
    I_LD1  = normalized_brightness_LD(rho, diam_LD1, u_lambda)
    
    I_tot = F_tot/(1+flux_ratio_LD1+ flux_ratio_G1)*I_UD*(1+flux_ratio_LD1*I_LD1/I_UD+flux_ratio_G1*I_G1/I_UD)  
    
    return rho, I_tot
