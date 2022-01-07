#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:03:17 2021

@author: jdrevon
"""
import numpy as np
from pretty_table_reading import READ_PRETTY_TABLE
from brightness_distribution_models import normalized_brightness_R, normalized_brightness_UD


def increase_intensity_profile_res(PATH_INTENSITY_PROFILE,**kwargs):

    
        header, data = READ_PRETTY_TABLE(PATH_INTENSITY_PROFILE,1)
        
        wavel = (header[1:]).astype(np.float)
        ext_radius = data[:,0]
        inner_radius = (data[:,0]-np.diff(data[:,0])[0])
        
        min_dist        = kwargs.get('min_dist', min(inner_radius))    
        max_dist        = kwargs.get('max_dist', max(ext_radius))    
        nb_points       = kwargs.get('nb_points', len(ext_radius))    


        I_tot_norm_HR = np.zeros((len(wavel),nb_points))
        rho   = np.linspace(min_dist,max_dist,nb_points)

        count=[]
        for t in range(len(ext_radius)):                
            cond = np.logical_and(rho<=ext_radius[t],rho>inner_radius[t])
            try:
                count = np.append(count,len(np.where(cond==True)[0]))            
            except :
                count = np.append(count,0)           

        if rho[0] == 0:
            count[0]=count[0]+1

        
        for i in range(len(wavel)):

            I_tot_norm_HR[i] = np.repeat(data[:,i+1],count.astype(np.int))

                            
        return rho, I_tot_norm_HR, wavel

    # plt.figure()
    # plt.imshow(profilegrid2,extent=[-R_image,R_image,-R_image,R_image])
    # plt.xlabel(r'$\alpha$ [mas]')
    # plt.ylabel(r'$\delta$ [mas]')
    # plt.colorbar(label=r'Intensity Ratio [I_tot/I_star]')
    # plt.title('Image of the model at %.2f µm'%wavel_UD) 



# plt.figure()
# plt.imshow(data_cube[0],extent=[-image_boundary,image_boundary,-image_boundary,image_boundary])
# plt.xlabel(r'$\alpha$ [mas]')
# plt.ylabel(r'$\delta$ [mas]')
# plt.colorbar(label=r'Intensity Ratio [I_tot/I_star]')
# plt.title('Image of the model at %.2f µm'%wavel)   
