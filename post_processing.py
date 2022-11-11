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
import numpy as np
import istarmap
import multiprocessing as mp
from multiprocessing import Pool
import tqdm
from rhapsody_init import nprocs

def image(i, PATH_OUTPUT_INT, PATH_INTENSITY_PROFILE, wavel_model, R_image, resolution):
    

    # if image_resolution == 'default':
    #     rho, I_tot_norm_HR, wavel = increase_intensity_profile_res(PATH_INTENSITY_PROFILE)
    # else :
    #     try:
    #         rho, I_tot_norm_HR, wavel = increase_intensity_profile_res(PATH_INTENSITY_PROFILE, nb_points=int(image_resolution/2))
    #     except:
    #         rho, I_tot_norm_HR, wavel = increase_intensity_profile_res_1W(PATH_INTENSITY_PROFILE, nb_points=int(image_resolution/2))
                
    try:
        rho, I_tot_norm_HR, wavel = increase_intensity_profile_res(PATH_INTENSITY_PROFILE, nb_points=int(resolution/2))
    except:
        rho, I_tot_norm_HR, wavel = increase_intensity_profile_res_1W(PATH_INTENSITY_PROFILE, nb_points=int(resolution/2))
        
    x_model_image, image_model = imaging_model(rho,I_tot_norm_HR[i],R_image)

    plot_intensity_model_image(x_model_image, image_model, wavel_model, R_image, PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT)

    
    x_model_profile = rho
    y_model_profile = I_tot_norm_HR[i]

    
    plot_intensity_model_profile(x_model_profile, y_model_profile, wavel_model, PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT, xlim_max = R_image)    
    
    
    return wavel_model, x_model_image, image_model, x_model_profile, y_model_profile



def post_processing(PATH_OUTPUT_FIT_RES, PATH_OUTPUT_INT, PATH_INTENSITY_PROFILE, list_wavel, image_rec_windows, resolution):

    print('STARTING PLOTTING INTENSITY RADIAL PROFILES AND EQUIVALENT 2D IMAGES')

    intensity = []
    distance  = []
    stock_images = []
    x_model_image = [] 
    wavel = []

    with mp.Pool(processes=nprocs) as pool:
        iterable = [(i, PATH_OUTPUT_INT, PATH_INTENSITY_PROFILE, list_wavel[i], image_rec_windows, resolution) for i in range(len(list_wavel))]
            
        for result in tqdm.tqdm(pool.istarmap(image, iterable),
                           total=len(iterable)):

            intensity.append(result[4])
            stock_images.append(result[2])
            wavel.append(result[0])
            
        distance = result[3]
        x_model_image = result[1]
    
    image_to_data_cube(stock_images, x_model_image , x_model_image, wavel, PATH_OUTPUT_FIT_RES, 'image_datacube')

    print('END PLOTTING FACED-ON INTENSITY RADIAL PROFILES AND EQUIVALENT 2D IMAGES')
    
    return np.array(distance), np.array(intensity), np.array(wavel)
    
