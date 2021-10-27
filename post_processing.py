# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 20:09:48 2021

@author: jdrevon
"""

from increase_intensity import increase_intensity_profile_res
from data_cube import image_to_data_cube
from imaging_model import imaging_model
from plot_intensity_rec_model import plot_intensity_model_image, plot_intensity_model_profile

def post_processing(PATH_OUTPUT_INT, wavel_model, R_image, image_resolution='default'):
    
    
    PATH_INTENSITY_PROFILE_LM = PATH_OUTPUT_INT+'intensity_LM.dat'
    PATH_INTENSITY_PROFILE_N  = PATH_OUTPUT_INT+'intensity_N.dat'

    if image_resolution == 'default':
        rho_LM, I_tot_norm_HR_LM, wavel_LM = increase_intensity_profile_res(PATH_INTENSITY_PROFILE_LM)
        rho_N, I_tot_norm_HR_N, wavel_N    = increase_intensity_profile_res(PATH_INTENSITY_PROFILE_N)
    else :
        rho_LM, I_tot_norm_HR_LM, wavel_LM = increase_intensity_profile_res(PATH_INTENSITY_PROFILE_LM, nb_points=int(image_resolution/2))
        rho_N, I_tot_norm_HR_N, wavel_N    = increase_intensity_profile_res(PATH_INTENSITY_PROFILE_N, nb_points=int(image_resolution/2))

    
    index_LM = 0
    index_N  = 0
    
    stock_images_LM = []
    stock_images_N  = []
    
    print(' Processing the %i wavelengths...'%(len(wavel_model)))

    for i in range(len(wavel_model)):
        
        
        if wavel_model[i]<6:
            if index_LM == 0:
                print(' Processing the LM-Band')
            
            x_model_image_LM, image_model_LM = imaging_model(rho_LM,I_tot_norm_HR_LM[index_LM],R_image)
            x_model_image, image_model = x_model_image_LM, image_model_LM     
            x_model_profile = rho_LM
            y_model_profile = I_tot_norm_HR_LM[index_LM]
            index_LM += 1        
            stock_images_LM.append(image_model_LM)
            
        else : 
            if index_N == 0:
                print(' Processing the N-Band')

            x_model_image_N, image_model_N = imaging_model(rho_N,I_tot_norm_HR_N[index_N],R_image)    
            x_model_image,image_model = x_model_image_N, image_model_N
            x_model_profile = rho_N
            y_model_profile = I_tot_norm_HR_N[index_N]
            index_N += 1        
            stock_images_N.append(image_model_N)


        # print('IMAGE START')
        plot_intensity_model_image(x_model_image, image_model, wavel_model[i], R_image, PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT)
        # print('IMAGE END')
        
        # print('PROFILE START')
        plot_intensity_model_profile(x_model_profile, y_model_profile, wavel_model[i], PLOT= False, SAVE_OUTPUT = PATH_OUTPUT_INT, xlim_max = R_image)
        # print('PROFILE END')

    image_to_data_cube(stock_images_LM, x_model_image_LM , x_model_image_LM, wavel_LM, PATH_OUTPUT_INT, 'image_LM_datacube')
    image_to_data_cube(stock_images_N, x_model_image_N , x_model_image_N, wavel_N, PATH_OUTPUT_INT, 'image_N_datacube')

    print('END IMAGE RECONSTRUCTION')
    
    
    
    return
