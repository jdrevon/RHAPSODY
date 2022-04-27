#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:14:23 2020

@author: Julien
"""

from rhapsody_init import OIFITS_FLUX, ERROR_SUP, DATA_DIR, DATA_band_name, PROCESS_DIR, REG_method, HP, FITTING, model_q, model_q_max, model_q_min
from initialisation_rings import rings
from fits_reading_dico import OIFITS_READING, OIFITS_SORTING
from stock_dico_values import stock_V2_from_dico
from create_folder import FOLDER_CREATION
from remove_files_from_folder import REMOVE_ALL_FILES_FROM_FOLDER
from fitting_function_ud_ring import UD_modeling
import numpy as np
from spectra_extraction import spectra_data
from data_copy import COPY_DATA
import sys
from big_spectra import big_spectra_with_flux,big_spectra_without_flux
from data_flag_flux import FLUX_FLAG
from model_visibilities import V_ring
from post_processing import post_processing


if __name__ == '__main__':
    
    # COPY THE DATA IN THE RIGHT FOLDER 

    COPY_DATA(DATA_DIR, PROCESS_DIR)

    # SORT DATA in FOLDERS based on V2 and T3 flags:
        
    OIFITS_SORTING()
        
    # READ ALL THE DATA and put them in dicos in order to manipulate the data easily
    
    OIFITS_TOT = OIFITS_READING()
        
    
    # We regroup all the usefull information that we want in simple masked array to handle easily the data for the fitting
    
    wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR = stock_V2_from_dico(OIFITS_TOT)

    # Let's correct the ERROR BARS here:

    for i in range(len(DATA_band_name)):
        
        V2_DATA_ERR[i] = np.sqrt(V2_DATA_ERR[i]**2+ERROR_SUP[i]**2)
        
    
    # Rings initialization
        
    diam_inner_ring, diam_outter_ring, I_norm, flux, flux_max = rings()
    
    # Pre-Computing for each wavelengths the associated modeled visibilities for each rings
    
    V_model = [] # I stock the visibilities for each band
    
    
    for i in range(len(DATA_band_name)):
        
        V_model_tmp=[]
        
        wavel_in_bdth = np.unique(wavel_DATA[i])
        
        for k in range(len(wavel_in_bdth)):
        
            cond = np.where(wavel_DATA[i]==wavel_in_bdth[k])
            q_model = np.array(q_DATA[i][cond])      
    
            V_tmp = [V_ring(q_model,diam_inner_ring[i][j],diam_outter_ring[i][j]) for j in range(len(diam_outter_ring[i]))]            
                
            V_model_tmp.append(V_tmp)
    
        V_model.append(V_model_tmp)

    # V_model strucutre : 
    # (N,M,K,L) : 
    # N = # of bandwidth, 
    # M = # of unique wavelength in the bandwidth
    # K = # of rings computed for the specific bandwidth 
    # L = # of visibilities computed for each ring at each spatial frequency
    
    
    
    # List of the spatial frequencies used to compute the model plots with higer resolution than the one provided by the instrument

    q_UD_HR = []
    
    for i in range(len(DATA_band_name)):
    
        baseline = q_DATA[i]*wavel_DATA[i]    
        baseline_max  = np.max(baseline)
        baseline_min  = np.min(baseline)

        if model_q_min[i]== None:
            model_q_min[i] = np.log10(baseline_min/np.max(wavel_DATA[i])) 
        if model_q_max[i]== None:
            model_q_max[i] = np.log10(baseline_max/np.min(wavel_DATA[i]))
        q_HR = np.logspace(model_q_min[i],model_q_max[i],model_q[i])
        q_UD_HR.append(q_HR)
    
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
    
        if FITTING == True:
            FOLDER_CREATION(PATH_RESULTS)
            FOLDER_CREATION(PATH_OUTPUT)
            REMOVE_ALL_FILES_FROM_FOLDER(PATH_OUTPUT)
            for k in range(len(DATA_band_name)):
                FOLDER_CREATION(PATH_OUTPUT_VIS[k]) 
                FOLDER_CREATION(PATH_OUTPUT_INT[k])
                FOLDER_CREATION(PATH_OUTPUT_HIST[k])
                FOLDER_CREATION(PATH_OUTPUT_FIT_RES[k])
                FOLDER_CREATION(PATH_OUTPUT_SPECTRA[k])
    
    
            # UD_modeling(wavel_ALL,wavel_MATISSE, q_MATISSE, V2_MATISSE, V2_MATISSE_ERR,\
            #                          q_UD_HR, diam_outter_ring, diam_inner_ring, flux_ratio_LM, flux_ratio_N,\
            #                                      HP[i], V_model,\
            #                                          PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)
    
            UD_modeling(wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR,\
                                     q_UD_HR, diam_outter_ring, diam_inner_ring, flux, I_norm, flux_max, \
                                                 HP[i], V_model,\
                                                     PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)
    
    
    

        # Plots: 
        
        for k in range(len(DATA_band_name)):
            
            list_wavel = np.unique(wavel_DATA[k])*1E6
            
            # 1/ Construction of the intensity profile plots + 2/ Construction of the image from the intensity profile plots
            
            PATH_INTENSITY_PROFILE = PATH_OUTPUT_FIT_RES[k]+'intensity_%s_band.dat'%DATA_band_name[k]

            
            distance, intensity, wavel = post_processing(PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_INT[k], PATH_INTENSITY_PROFILE, list_wavel, 40, image_resolution=2**8)    
    
            # RESCALE 
            
            W,D=np.meshgrid(np.append(wavel,wavel[-1]+np.diff(wavel)[-1]),np.append([0],distance))
    
            # PLOT BIG SPECTRA:
                
            # WITH FLUX :

            if OIFITS_FLUX == True:
                # PROCESSING AND FLAGGING THE FLUX    
                print('STARTING FLUX SORTING')
                # Here I flag the flu data 
                OIFITS_TOT = FLUX_FLAG(OIFITS_TOT)
                # Here I compute the spectra for each AT's and the total spectra
                FLUX_WAVEL, FLUX_DATA, FLUX_DATA_err = spectra_data(OIFITS_TOT[k], PATH_OUTPUT_FIT_RES[k], PATH_OUTPUT_SPECTRA[k])
                # Here I plot the BIG SPECTRA
                big_spectra_with_flux(W,D,intensity.T, FLUX_WAVEL, FLUX_DATA, PATH_OUTPUT_SPECTRA[k])
                
            # WITHOUT FLUX:
    
            else:
                big_spectra_without_flux(W,D,intensity.T, PATH_OUTPUT_SPECTRA[k])

