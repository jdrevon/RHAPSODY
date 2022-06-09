#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 20:09:48 2021

@author: jdrevon
"""

from increase_intensity import increase_intensity_profile_res,increase_intensity_profile_res_1W
from data_cube import image_to_data_cube
from imaging_model import imaging_model
from plot_intensity_rec_model import plot_intensity_model_image, plot_intensity_model_profile

def post_processing(PATH_OUTPUT_FIT_RES, PATH_OUTPUT_INT, PATH_INTENSITY_PROFILE, wavel_model, R_image, image_resolution='default'):
    
    print('STARTING PLOTTING INTENSITY RADIAL PROFILES AND EQUIVALENT 2D IMAGES')

    if image_resolution == 'default':
        rho, I_tot_norm_HR, wavel = increase_intensity_profile_res(PATH_INTENSITY_PROFILE)
    else :
        try:
            rho, I_tot_norm_HR, wavel = increase_intensity_profile_res(PATH_INTENSITY_PROFILE, nb_points=int(image_resolution/2))
        except:
            rho, I_tot_norm_HR, wavel = increase_intensity_profile_res_1W(PATH_INTENSITY_PROFILE, nb_points=int(image_resolution/2))
                
    stock_images = []    

    for i in range(len(wavel_model)):        

        print('Wavelength n°%i/%i'%(i+1,len(wavel_model)))
        
        x_model_image, image_model = imaging_model(rho,I_tot_norm_HR[i],R_image)
        stock_images.append(image_model)

        plot_intensity_model_image(x_model_image, image_model, wavel_model[i], R_image, PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT)

        
        x_model_profile = rho
        y_model_profile = I_tot_norm_HR[i]

        
        plot_intensity_model_profile(x_model_profile, y_model_profile, wavel_model[i], PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT, xlim_max = R_image)

    image_to_data_cube(stock_images, x_model_image , x_model_image, wavel, PATH_OUTPUT_FIT_RES, 'image_datacube')

    print('END PLOTTING INTENSITY RADIAL PROFILES AND EQUIVALENT 2D IMAGES')
    
    
    
    return rho, I_tot_norm_HR, wavel
