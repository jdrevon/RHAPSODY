#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:14:23 2020

@author: Julien
"""

from rhapsody_init import OIFITS_FLUX, ERROR_SUP, DATA_DIR, DATA_band_name, PROCESS_DIR, REG_method, HP, FITTING, model_visibilities, image_rec_windows, inc_flag, READING_ONLY, borders_spectra, manual_borders, additionnal_ticks, CONCATENATE, resolution, fft_q_min, fft_q_max, fft_nb_points
from initialisation_rings import rings
from fits_reading_dico import OIFITS_READING, OIFITS_SORTING, OIFITS_READING_concatenate
from stock_dico_values import stock_V2_from_dico, stock_V2_from_dico_CONCATENATE
from create_folder import FOLDER_CREATION
from remove_files_from_folder import REMOVE_ALL_FILES_FROM_FOLDER
from fitting_function_ud_ring import UD_modeling
import numpy as np
from spectra_extraction import spectra_data, spectra_data_CONCATENATE
from data_copy import COPY_DATA
import sys
from big_spectra import big_spectra_with_flux,big_spectra_without_flux
from data_flag_flux import FLUX_FLAG, FLUX_FLAG_CONCATENATE
from model_visibilities import V_ring
from post_processing import post_processing
from numpy import cos as cos
from image_reconstruction import image_reconstruction

if __name__ == '__main__':
    
    
    if READING_ONLY == False:
        
        user = input("WARNING! Previously computed data with the same band name located in the following directory will be erased:\n%s.\nDo You Still Want To Continue? [y/n] \n"%PROCESS_DIR)
        
        if user == "n":
            print('Wise decision... RHAPSODY has shut down')
            sys.exit()
        elif user == "y":
            print("OK let's goooooo!")
        
        else:
            print("Mamaaaaaa Ooooooooo !!!!! Let's go listen Bohemian Rhapsody and come back later (next time, please chose y or n as an answer please)")
            sys.exit()
        
    # COPY THE DATA IN THE RIGHT FOLDER 

    COPY_DATA(DATA_DIR, PROCESS_DIR)

    # SORT DATA in FOLDERS based on V2 and T3 flags:
    
    if CONCATENATE == False:
        
        OIFITS_SORTING()
        
        # READ ALL THE DATA and put them in dicos in order to manipulate the data easily
        OIFITS_TOT = OIFITS_READING()
        wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR, U_DATA, V_DATA = stock_V2_from_dico(OIFITS_TOT)
    
    elif CONCATENATE == True:
        
        OIFITS_SORTING()
        OIFITS_TOT = OIFITS_READING_concatenate()
        wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR, U_DATA, V_DATA = stock_V2_from_dico_CONCATENATE(OIFITS_TOT)
        
    
    # We regroup all the usefull information that we want in simple masked array to handle easily the data for the fitting
    

    # Let's correct the ERROR BARS here:

    for i in range(len(DATA_band_name)):
        
        V2_DATA_ERR[i] = np.sqrt(V2_DATA_ERR[i]**2+ERROR_SUP[i]**2)
        
    qu_DATA = U_DATA/wavel_DATA
    qv_DATA = V_DATA/wavel_DATA
    
    # Rings initialization
        
    diam_inner_ring, diam_outter_ring, I_norm, flux, flux_max = rings()
    
    # Pre-Computing for each bandwidth the associated modeled visibilities for each rings
    
    V_model = [] # I stock the visibilities for each band
    qu_interp = []
    qv_interp = []
    q_interp  = []


    for i in range(len(DATA_band_name)):
        
                
        qu_DATA_max = np.log10(np.amax(np.abs(qu_DATA[i]))) # +0.1 is here to avoid interpolation problem
        qv_DATA_max = np.log10(np.amax(np.abs(qv_DATA[i])))#*1/cos(np.deg2rad(87)))
        
        qu_interp_tmp = np.array([0]+np.logspace(4, qu_DATA_max, model_visibilities[i]).tolist())
        qv_interp_tmp = np.array([0]+np.logspace(4, qv_DATA_max, model_visibilities[i]).tolist())

        V_tmp = [V_ring(qu_interp_tmp, qv_interp_tmp, diam_inner_ring[i][j],diam_outter_ring[i][j]) for j in range(len(diam_outter_ring[i]))]            
                    
        qu_interp.append(qu_interp_tmp)
        qv_interp.append(qv_interp_tmp)     
        q_interp.append(np.sqrt(qu_interp_tmp**2+qv_interp_tmp**2))
        
        V_model.append(V_tmp)

    # V_model strucutre : 
    # (N,K,L) : 
    # N = # of bandwidth, 
    # K = # of rings computed for the specific bandwidth 
    # L = # of visibilities computed for each ring at each spatial frequency
    
    
    # Initialization of the RESULTS folder to put model results
    
    PATH_RESULTS = PROCESS_DIR+'/RESULTS'+ '/'
    
    if REG_method == 'TV':    
    
        NAME_FOLDER_tmp = 'REG_TV_HP_'
    
    elif REG_method == 'QS':
    
        NAME_FOLDER_tmp = 'REG_TV_HP_'
    
    else:
        print("Please enter a valid regularization in the file 'RHAPSODY_init' : 'TV' for Total Variation or 'QS' for Quadratic Smoothness")    
        sys.exit()    

    HP_string = ['%.1E' %HP[i] for i in range(len(HP))]
    
    for i in range(len(HP)):
        for k in range(len(DATA_band_name)):
            PATH_OUTPUT = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] 
            PATH_OUTPUT_VIS = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/VISIBILITY/' for k in range(len(DATA_band_name))]
            PATH_OUTPUT_INT = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/INTENSITY/' for k in range(len(DATA_band_name))]
            PATH_OUTPUT_HIST = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/HISTOGRAM/' for k in range(len(DATA_band_name))]
            PATH_OUTPUT_FIT_RES = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/FIT_RESULTS/' for k in range(len(DATA_band_name))]
            PATH_OUTPUT_SPECTRA = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/SPECTRA/' for k in range(len(DATA_band_name))]
            if inc_flag==True:
                PATH_OUTPUT_IMAGE_REC = [PATH_OUTPUT + '/' + DATA_band_name[k] + '/IMAGE_REC/' for k in range(len(DATA_band_name))]
    
        if READING_ONLY == False:            
            
            FOLDER_CREATION(PATH_RESULTS)
            FOLDER_CREATION(PATH_OUTPUT)
            
            #If a folder with the same band name is still here, delete it (DO IT!)
            try:
                for k in range(len(DATA_band_name)):
                    REMOVE_ALL_FILES_FROM_FOLDER(PATH_OUTPUT+ '/' +DATA_band_name[k])
            except:
                None
            for k in range(len(DATA_band_name)):
                FOLDER_CREATION(PATH_OUTPUT_VIS[k]) 
                FOLDER_CREATION(PATH_OUTPUT_INT[k])
                FOLDER_CREATION(PATH_OUTPUT_HIST[k])
                FOLDER_CREATION(PATH_OUTPUT_FIT_RES[k])
                FOLDER_CREATION(PATH_OUTPUT_SPECTRA[k])
                if inc_flag==True:
                    FOLDER_CREATION(PATH_OUTPUT_IMAGE_REC[k])
    
    
            # UD_modeling(wavel_ALL,wavel_MATISSE, q_MATISSE, V2_MATISSE, V2_MATISSE_ERR,\
            #                          q_UD_HR, diam_outter_ring, diam_inner_ring, flux_ratio_LM, flux_ratio_N,\
            #                                      HP[i], V_model,\
            #                                          PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)
    
            UD_modeling(wavel_DATA, q_DATA, qu_DATA, qv_DATA, V2_DATA, V2_DATA_ERR,\
                                     q_interp, qu_interp, qv_interp, diam_outter_ring, diam_inner_ring, flux, I_norm, flux_max, \
                                                 HP[i], V_model,\
                                                     PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)
    
    
    

        # Plots: 
        
        for k in range(len(DATA_band_name)):
            
            list_wavel = np.unique(wavel_DATA[k])*1E6
            
            # 1/ Construction of the intensity profile plots + 2/ Construction of the image from the intensity profile plots
            
            PATH_INTENSITY_PROFILE = PATH_OUTPUT_FIT_RES[k]+'intensity_%s_band.dat'%DATA_band_name[k]

            
            distance, intensity, wavel = post_processing(PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_INT[k], PATH_INTENSITY_PROFILE, list_wavel, image_rec_windows[k], resolution[k])    
    
            # SPECTRA
        
            # Apply the spectra only if we have more than 1 WVL to model
            
            if np.shape(list_wavel)[0] != 1:
                W,D=np.meshgrid(np.append(wavel,wavel[-1]+np.diff(wavel)[-1]),np.append(distance,max(distance)+np.diff(distance)[0]))

                # WITH FLUX :
    
                if OIFITS_FLUX == True:
                    # PROCESSING AND FLAGGING THE FLUX    
                    print('STARTING FLUX SORTING')
                    # Here I flag the flu data 
                    if CONCATENATE == True :
                        OIFITS_TOT = FLUX_FLAG_CONCATENATE(OIFITS_TOT)
                        FLUX_WAVEL, FLUX_DATA, FLUX_DATA_err = spectra_data_CONCATENATE(OIFITS_TOT[k], PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_SPECTRA[k])
                        
                    else:                            
                        OIFITS_TOT = FLUX_FLAG(OIFITS_TOT)
                        FLUX_WAVEL, FLUX_DATA, FLUX_DATA_err = spectra_data(OIFITS_TOT[k], PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_SPECTRA[k])
                    # Here I compute the spectra for each AT's and the total spectra
                    

                    # Here I plot the BIG SPECTRA
                    if manual_borders == True:
                        big_spectra_with_flux(W,D,intensity.T, FLUX_WAVEL, FLUX_DATA, PATH_OUTPUT_SPECTRA[k], additionnal_ticks[k], borders = borders_spectra[k])
                    else:
                        big_spectra_with_flux(W,D,intensity.T, FLUX_WAVEL, FLUX_DATA, PATH_OUTPUT_SPECTRA[k], additionnal_ticks[k])

                # WITHOUT FLUX:
        
                else:
                    if manual_borders == True:
                        big_spectra_without_flux(W,D,intensity.T, PATH_OUTPUT_SPECTRA[k], additionnal_ticks[k],  borders = borders_spectra[k])
                    else: 
                        big_spectra_without_flux(W,D,intensity.T, PATH_OUTPUT_SPECTRA[k], additionnal_ticks[k])

            #IMAGE RECONSTRUCTION:
            if inc_flag == True:
                image_reconstruction(PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_IMAGE_REC[k], DATA_band_name[k], diam_inner_ring[k], diam_outter_ring[k], image_rec_windows[k], resolution[k], fft_q_min[k], fft_q_max[k], fft_nb_points[k])
            else: None
                

                
