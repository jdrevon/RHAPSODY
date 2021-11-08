# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:14:23 2020

@author: Julien
"""

from rhapsody_init import size_ring_LM, size_ring_N, NBR_ring_LM, NBR_ring_N, alpha_LM, alpha_N, ERROR_SUP, MATISSE_DIR, REG_method, HP, FITTING, model_q_LM, model_q_N, model_q_max_LM, model_q_max_N, model_q_min_LM, model_q_min_N, BB_norm, BB_temperature, stellar_radii, distance_target
from initialisation_rings import rings
from fits_reading_dico import OIFITS_READING, OIFITS_SORTING
from stock_dico_values import stock_V2_from_dico
from create_folder import FOLDER_CREATION
from remove_files_from_folder import REMOVE_ALL_FILES_FROM_FOLDER
from fitting_function_ud_ring import UD_modeling
import numpy as np
from spectra_extraction import spectra_matisse
from A1_mas_to_rad import au_to_R_sun
from black_body import BB_v_Jy


import sys

from model_visibilities import V_ring

from output_intensity_profile import output_intensity_profile

from post_processing import post_processing

# SORT DATA in FOLDERS based on V2 and T3 flags:
    
OIFITS_SORTING()

# READ ALL THE DATA and put them in two dicos in order to manipulate the data easily

OIFITS_TOT_LM, OIFITS_TOT_N = OIFITS_READING()

# We regroup all the usefull information that we want in simple masked array to handle easily the data for the fitting

wavel_MATISSE, q_MATISSE, V2_MATISSE, V2_MATISSE_ERR = stock_V2_from_dico(OIFITS_TOT_LM, OIFITS_TOT_N)

wavel_LM = np.unique(wavel_MATISSE[wavel_MATISSE<6*1E-6])*1E6 #µm
wavel_N  = np.unique(wavel_MATISSE[wavel_MATISSE>6*1E-6])*1E6 #µm

wavel_ALL = np.concatenate((wavel_LM,wavel_N)) #µm

# Rings initialization


diam_inner_ring, diam_outter_ring, I_norm_LM, I_norm_N, flux_LM, flux_N, flux_LM_max, flux_N_max = rings(size_ring_LM, size_ring_N, NBR_ring_LM, NBR_ring_N, len(wavel_LM), len(wavel_N), alpha_LM, alpha_N)

# Pre-Computing for each wavelengths the associated modeled visibilities for each rings

V_model = []

for i in range(len(wavel_ALL)):

    cond = np.where(wavel_MATISSE*1E6==wavel_ALL[i])
    q_model = np.array(q_MATISSE[cond])      

    if wavel_ALL[i]<6:    
        V_tmp = [V_ring(q_model,diam_inner_ring[0][k],diam_outter_ring[0][k]) for k in range(len(diam_outter_ring[0]))]            
        
    else : 
        V_tmp = [V_ring(q_model,diam_inner_ring[1][k],diam_outter_ring[1][k]) for k in range(len(diam_outter_ring[1]))]

    V_model.append(V_tmp)



# List of the spatial frequency used to compute the model plots with higer resolution than the one provided by the instrument

baseline_LM = q_MATISSE[wavel_MATISSE<6E-6]*wavel_MATISSE[wavel_MATISSE<6E-6]
baseline_N  = q_MATISSE[wavel_MATISSE>6E-6]*wavel_MATISSE[wavel_MATISSE>6E-6]

baseline_max_LM  = np.max(baseline_LM)
baseline_min_LM  = np.min(baseline_LM)

baseline_max_N   = np.max(baseline_N)
baseline_min_N   = np.min(baseline_N)


if model_q_min_LM== None:
    model_q_min_LM = np.log10(baseline_min_LM/np.max(wavel_MATISSE[wavel_MATISSE<6E-6])) 

if model_q_max_LM== None:
    model_q_max_LM = np.log10(baseline_max_LM/np.min(wavel_MATISSE[wavel_MATISSE<6E-6]))

if model_q_min_N == None:
    model_q_min_N = np.log10(baseline_min_N/np.max(wavel_MATISSE[wavel_MATISSE>6E-6])) 

if model_q_max_N == None:
    model_q_max_N = np.log10(baseline_max_N/np.min(wavel_MATISSE[wavel_MATISSE>6E-6]))


q_UD_HR = np.array([np.logspace(model_q_min_LM,model_q_max_LM,model_q_LM), np.logspace(model_q_min_N,model_q_max_N,model_q_N)])

    
V2_MATISSE_ERR = np.sqrt(V2_MATISSE_ERR**2+ERROR_SUP**2)


# Initialization of the RESULTS folder to put model results

PATH_RESULTS = MATISSE_DIR+'/RESULTS'+ '/'

if REG_method == 'TV':    

    NAME_FOLDER_tmp = 'REG_TV_HP_'

elif REG_method == 'QS':

    NAME_FOLDER_tmp = 'REG_TV_HP_'

else:
    print("Please enter a valid regularization in the file 'RHAPSODY_init' : 'TV' for Total Variation or 'QS' for Quadratic Smoothness")    
    sys.exit()    

HP_string = ['%.1E' %HP[i] for i in range(len(HP))]

for i in range(len(HP)):
    PATH_OUTPUT = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i]
    PATH_OUTPUT_VIS = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'VISIBILITY/'
    PATH_OUTPUT_INT = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'INTENSITY/'
    PATH_OUTPUT_HIST = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'HISTOGRAM/'
    PATH_OUTPUT_FIT_RES = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'FIT_RESULTS/'
    PATH_OUTPUT_SPECTRA = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'SPECTRA/'

    if FITTING == True:
        FOLDER_CREATION(PATH_RESULTS)
        FOLDER_CREATION(PATH_OUTPUT)
        REMOVE_ALL_FILES_FROM_FOLDER(PATH_OUTPUT)
        FOLDER_CREATION(PATH_OUTPUT_VIS)
        FOLDER_CREATION(PATH_OUTPUT_INT)
        FOLDER_CREATION(PATH_OUTPUT_HIST)
        FOLDER_CREATION(PATH_OUTPUT_FIT_RES)
        FOLDER_CREATION(PATH_OUTPUT_SPECTRA)


        # UD_modeling(wavel_ALL,wavel_MATISSE, q_MATISSE, V2_MATISSE, V2_MATISSE_ERR,\
        #                          q_UD_HR, diam_outter_ring, diam_inner_ring, flux_ratio_LM, flux_ratio_N,\
        #                                      HP[i], V_model,\
        #                                          PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)

        UD_modeling(wavel_ALL,wavel_MATISSE, q_MATISSE, V2_MATISSE, V2_MATISSE_ERR,\
                                 q_UD_HR, diam_outter_ring, diam_inner_ring, flux_LM, flux_N, I_norm_LM, I_norm_N, flux_LM_max, flux_N_max, \
                                             HP[i], V_model,\
                                                 PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES,  PLOT=False)


        # Plots: 
            
        # 1/ Construction of the intensity profile plots + 2/ Construction of the image from the intensity profile plots
        
        post_processing(PATH_OUTPUT_INT, wavel_ALL, 40, image_resolution=2**8)


        # RESCALE 
        
        distance_LM, intensity_LM, distance_N, intensity_N, wavel_LM, wavel_N = output_intensity_profile(PATH_OUTPUT_INT+'intensity_LM.dat', PATH_OUTPUT_INT+'intensity_N.dat',wavel_ALL)
        W_LM,D_LM=np.meshgrid(np.append(wavel_LM,wavel_LM[-1]+np.diff(wavel_LM)[-1]),np.append([0],distance_LM))
        W_N,D_N=np.meshgrid(np.append(wavel_N,wavel_N[-1]+np.diff(wavel_N)[-1]),np.append([0],distance_N))
    
        # RESCALE 
        widths  = [1,1]
        heights = [1,5]
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        import matplotlib
        import matplotlib.colors as colors
        from cmcrameri import cm
    
        mp   = cm.lajolla
        rmap = mp.reversed()     
    
        fig=plt.figure(figsize=(8, 12))
        gs1 = GridSpec(2, 2, width_ratios=widths, height_ratios=heights)
        
    
        ax1 = plt.subplot(gs1[1,0])
        ax2 = plt.subplot(gs1[1,1])
    
        ax3 = plt.subplot(gs1[0,0])
        ax4 = plt.subplot(gs1[0,1])
    
        
        plot1=ax1.pcolor(W_LM,D_LM,np.array(intensity_LM),norm=colors.LogNorm(vmin=1E-4, vmax=1),shading='auto', cmap=rmap)
        
        W_LM_2,D_LM_2=np.meshgrid(wavel_LM,distance_LM)
    
        ax1.minorticks_on()
        rect=plt.Rectangle((3.05, 0), 0.25, 56,alpha=0.4, color='lightgrey')
        ax1.add_patch(rect)
        rect2=plt.Rectangle((3.75, 0), 0.25, 56,alpha=0.4, color='lightgrey')
        ax1.add_patch(rect2)
        ax1.set_yscale('log')
        ax1.set_ylim(3,max(distance_N))
        ax1.yaxis.set_ticks_position('both')
    
        ax1.set_ylabel('Angular distance from the star center [mas]',fontsize = 13)
        ax1.set_yticks([5,10,12,15,18,25,30,40,50,75,100])    
        ax1.tick_params(axis='both', labelsize=13)
        ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
        ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
        plot2=ax2.pcolor(W_N,D_N,np.array(intensity_N),norm=colors.LogNorm(vmin=1E-4, vmax=1), cmap=rmap)
        
        W_N_2,D_N_2=np.meshgrid(wavel_N,distance_N)
    
        ax2.set_ylim(3,max(distance_N))
        ax2.minorticks_on()
        ax2.set_yscale('log')
        ax2.set_yticks([5,10,12,15,18,25,30,40,50,75,100])    
        ax2.yaxis.set_ticks_position('both')
        ax2.tick_params(axis='both', labelsize=13)
        ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.96, 0.15, 0.05, 0.7])
        fig.colorbar(plot1, cax=cbar_ax).set_label(label=r'$I(\rho)/I(0)$',size=13)
        fig.colorbar(plot1, cax=cbar_ax).ax.tick_params(labelsize=13) 

        
        FLUX_WAVEL, FLUX_DATA, FLUX_DATA_err = spectra_matisse(OIFITS_TOT_LM, OIFITS_TOT_N)

        wavel_MATISSE_flux = np.array(FLUX_WAVEL)    
        

        if BB_norm == True:
        
            wavel_BB, blackbody   = BB_v_Jy(wavel_MATISSE_flux*1E-6, BB_temperature, au_to_R_sun(stellar_radii), distance_target)
            data2                 = FLUX_DATA/blackbody.astype('float')
            MATISSE_flux_LM       = data2/np.mean(data2[np.logical_and(wavel_MATISSE_flux>2,wavel_MATISSE_flux<5)])
            MATISSE_flux_N        = data2/np.mean(data2[np.logical_and(wavel_MATISSE_flux>8,wavel_MATISSE_flux<12)])

        else:

            data2                 = FLUX_DATA            
            MATISSE_flux_LM       = data2/np.mean(data2[np.logical_and(wavel_MATISSE_flux>2,wavel_MATISSE_flux<5)])
            MATISSE_flux_N        = data2/np.mean(data2[np.logical_and(wavel_MATISSE_flux>8,wavel_MATISSE_flux<12)])

        
        ax3.tick_params(axis='both', labelsize=15)
        ax3.set_xlim(np.amin(W_LM),np.amax(W_LM))
        ax3.set_ylim(0.5,1.5)
        ax3.set_ylabel('Normalized Flux',fontsize =15)
    
        ax3.plot(np.sort(wavel_MATISSE_flux),[x for _, x in sorted(zip(wavel_MATISSE_flux, MATISSE_flux_LM))], label = 'MATISSE')
        ax3.minorticks_on()
        ax3.set
    
        ax4.tick_params(axis='both', labelsize=15)
        ax4.set_xlim(np.amin(W_N),np.amax(W_N))
        ax4.set_ylim(0.5,1.5)
        ax4.plot(np.sort(wavel_MATISSE_flux),[x for _, x in sorted(zip(wavel_MATISSE_flux, MATISSE_flux_N))])
        ax4.minorticks_on()

    
        handles, labels = [(a + b + c) for a, b, c in zip(ax1.get_legend_handles_labels(), ax2.get_legend_handles_labels(), ax3.get_legend_handles_labels())]
        fig1=ax1.figure
        fig1.text(0.5,0.04, 'Wavelength [µm]',fontsize = 15, ha="center", va="center")
        
        plt.setp(ax2.get_yticklabels(), visible=False)
    
        plt.savefig(PATH_OUTPUT_SPECTRA+'big_spectra.jpg',bbox_inches='tight', dpi=400)
        plt.close(fig)
