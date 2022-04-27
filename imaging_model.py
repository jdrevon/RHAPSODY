#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:39:37 2021

@author: jdrevon
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate


def imaging_model(rho,I_profile, R_image, **kwargs):
    
    method        = kwargs.get('method', 'linear')    

    
    xv = np.linspace(-R_image, R_image, len(rho)*2, endpoint=False) # Taille de la fenetre

    interpol_index = interp1d(rho, I_profile,kind=method)    
    X, Y = np.meshgrid(xv, xv)
    profilegrid2 = np.zeros(X.shape, float)
    current_radius = np.sqrt(X**2 + Y**2)
    cond=np.logical_and(current_radius<=max(rho),current_radius>=min(rho)) # Min et max des données
    profilegrid2[cond] = interpol_index(current_radius[cond])    

    return xv, profilegrid2  


def image_increase_resolution(x,y,image,factor):
        
    XX,YY = np.meshgrid(x,y)
    bound_x = min(abs(min(x)),max(x))
    bound_y = min(abs(min(y)),max(y))
    x_new = np.linspace(-bound_x,bound_y,len(x)*factor,endpoint=False)
    y_new = np.linspace(-bound_y,bound_y,len(y)*factor,endpoint=False)
    
    XX_new, YY_new = np.meshgrid(x_new,y_new)
    
    image_new = np.reshape(interpolate.griddata((XX.ravel(), YY.ravel()), image.ravel(), (XX_new.ravel(), YY_new.ravel()), method='linear'),np.shape(YY_new))

    # plt.figure()
    # plt.imshow(image_new,extent=[-bound_x,bound_x,-bound_y,bound_y])
    
    return x_new,y_new,image_new



def radial_profile_image(x,y,image,**kwargs):

    
    r_interp_max        = kwargs.get('r_interp_max', 40)    
    r_interp_min        = kwargs.get('r_interp_min', 0)    
    r_interp_nb         = kwargs.get('r_interp_nb', 500)    

    x0 = x[np.unique(np.where(image==np.amax(image))[0])[int(len(np.unique(np.where(image==np.amax(image))[0]))/2)]]
    y0 = y[np.unique(np.where(image==np.amax(image))[1])[int(len(np.unique(np.where(image==np.amax(image))[1]))/2)]]

    XX,YY = np.meshgrid((x-x0),(y-y0))
    
    TT_0 = np.arctan2(YY,XX)
    TT = np.rad2deg(TT_0)
    RR = np.sqrt(XX**2+YY**2)
    
    r_interp =  np.linspace(r_interp_min,r_interp_max,r_interp_nb,endpoint=False)
    theta_interp = np.linspace(-180,180,len(r_interp),endpoint=False) 
    
    RR_fit,TT_fit = np.meshgrid(r_interp,theta_interp)        

    interp = np.reshape(interpolate.griddata((RR.ravel(), TT.ravel()), image.ravel(), (RR_fit.ravel(), TT_fit.ravel()), method='linear'),np.shape(RR_fit))
    
    radial_profile = np.nanmean(interp,axis=0)
    std_deviation  = np.nanstd(interp,axis=0)        
   
    return r_interp, radial_profile, std_deviation

def radial_properties_image(x,y,image,radial_list):

    
    x0 = x[np.unique(np.where(image==np.amax(image))[0])[int(len(np.unique(np.where(image==np.amax(image))[0]))/2)]]
    y0 = y[np.unique(np.where(image==np.amax(image))[1])[int(len(np.unique(np.where(image==np.amax(image))[1]))/2)]]

    XX,YY = np.meshgrid((x-x0),(y-y0))
    
    TT_0 = np.arctan2(YY,XX)
    TT = np.rad2deg(TT_0)
    RR = np.sqrt(XX**2+YY**2)/1E7
    # r_interp =  np.linspace(-30,30,100)
    
    r_interp =  radial_list/1E7
    theta_interp = np.linspace(-180,180,len(r_interp),endpoint=False) 
    
    RR_fit,TT_fit = np.meshgrid(r_interp,theta_interp)        
    
    interp = np.reshape(interpolate.griddata((RR.ravel(), TT.ravel()), image.ravel(), (RR_fit.ravel(), TT_fit.ravel()), method='linear'),np.shape(RR_fit))
    
    radial_profile = np.nanmean(interp,axis=0)
    max_radial  = np.nanmax(interp,axis=0) - radial_profile  
    min_radial  = radial_profile-np.nanmin(interp,axis=0)  
    
   
    return r_interp*1E7, radial_profile, min_radial,max_radial


from brightness_distribution_models import normalized_brightness_UD, normalized_brightness_R
from pretty_table_reading import READ_PRETTY_TABLE
from A6_nearest_index import nearest_index

def increase_intensity_resolution(FILE_PATH_FLUX, len_UD, wavel_um):
    
        header,data = READ_PRETTY_TABLE(FILE_PATH_FLUX,1)

        
        rho_UD = header[1:-1].astype(np.float)
        res_UD = rho_UD*2
        res_UD_bef = res_UD-np.diff(res_UD)[0]


        nb_UD = len(rho_UD)

        rho_UD_new = np.linspace(0, max(rho_UD), len_UD)
        
        index = nearest_index(data[:,0], wavel_um)
        
        res_F_mod = data[index,1:]
        
        I_UD_norm = normalized_brightness_UD(rho_UD_new,res_UD[0])/normalized_brightness_UD(rho_UD_new,res_UD[0])[0]
        F_UD_0 = res_F_mod[0]
        I_tot_norm_OUT = I_UD_norm + np.sum([res_F_mod[k]/F_UD_0\
                                         *normalized_brightness_R(rho_UD_new, res_UD_bef[k], res_UD[k])/normalized_brightness_UD(rho_UD_new,res_UD[0])[0]\
                                         for k in range(1,nb_UD)],axis=0)

        return rho_UD_new, I_tot_norm_OUT    
    