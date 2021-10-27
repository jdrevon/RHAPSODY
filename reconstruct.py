# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:14:23 2020

@author: Julien
"""

from rhapsody_init import size_ring_LM, size_ring_N, NBR_ring_LM, NBR_ring_N, alpha_LM, alpha_N, ERROR_VIS, MATISSE_DIR, REG_method, HP, FITTING, model_q_LM, model_q_N, model_q_max_LM, model_q_max_N, model_q_min_LM, model_q_min_N
from initialisation_rings import rings
from fits_reading_dico import OIFITS_READING, OIFITS_SORTING
from stock_dico_values import stock_V2_from_dico
from create_folder import FOLDER_CREATION
from remove_files_from_folder import REMOVE_ALL_FILES_FROM_FOLDER
from fitting_function_ud_ring import UD_modeling

import numpy as np
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

diam_inner_ring, diam_outter_ring, I_norm_LM, I_norm_N, flux_ratio_LM, flux_ratio_N = rings(size_ring_LM, size_ring_N, NBR_ring_LM, NBR_ring_N, len(wavel_LM), len(wavel_N), alpha_LM, alpha_N)

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

if ERROR_VIS == True:
    
    V2_MATISSE_ERR = np.sqrt(V2_MATISSE_ERR**2+0.02**2)


# Initialization of the RESULTS folder to put model results

PATH_RESULTS = MATISSE_DIR+'/RESULTS'+ '/'
FOLDER_CREATION(PATH_RESULTS)

if REG_method == 'TV':    

    NAME_FOLDER_tmp = 'REG_TV_HP_'

elif REG_method == 'QS':

    NAME_FOLDER_tmp = 'REG_TV_HP_'

else:
    print("Please enter a valid regularization in the file 'RHAPSODY_init' : 'TV' for Total Variation or 'QS' for Quadratic Smoothness")    
    sys.exit()    

HP_string = ['%.1E' %HP[i] for i in range(len(HP))]

# test = [0.057978,0.146275,0.19457,0.200163,0.166722, 0.107445,0.04434,0.00234,1e-06,0.0, 0.000651,0.007195,0.010941,0.00791,0.001993,0.000568,0.00147, 4.2e-05, 0.006655, 0.01443, 0.00214, 1.6e-05, 0.0, 0.0, 0.0, 1e-06, 0.0, 0.0, 0.0,0.0,1e-06,1.2e-05,0.001483,0.00985,0.002095,0.000175, 1.5e-05,1e-06,2e-06,0.0,0.0,1e-06, 1e-06,0.0,0.0,0.0,0.0,0.0,0.0   ,0.0,1e-06,1e-06,0.000856,0.010257,0.001338, 3.6e-05]
# test = np.array(test)/max(test)

for i in range(len(HP)):
    PATH_OUTPUT = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i]
    PATH_OUTPUT_VIS = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'VISIBILITY/'
    PATH_OUTPUT_INT = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'INTENSITY/'
    PATH_OUTPUT_HIST = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'HISTOGRAM/'
    PATH_OUTPUT_FIT_RES = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'FIT_RESULTS/'
    PATH_OUTPUT_SPECTRA = PATH_RESULTS+NAME_FOLDER_tmp+HP_string[i] + '/' + 'SPECTRA/'

    if FITTING == True:
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
                                 q_UD_HR, diam_outter_ring, diam_inner_ring, flux_ratio_LM, flux_ratio_N,\
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
    
        handles, labels = [(a + b + c) for a, b, c in zip(ax1.get_legend_handles_labels(), ax2.get_legend_handles_labels(), ax3.get_legend_handles_labels())]
        fig1=ax1.figure
        fig1.text(0.5,0.04, 'Wavelength [µm]',fontsize = 15, ha="center", va="center")
        
        plt.setp(ax2.get_yticklabels(), visible=False)
    
        plt.savefig(PATH_OUTPUT_SPECTRA+'BIG_SPECTRA.jpeg',bbox_inches='tight', dpi=400)
        plt.close(fig)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        # RAPPORTS DE FLUX:
        
        # flux_ratio_PP_work(PATH_OUTPUT_tmp, NAME_FOLDER, OIFITS_TOT_LM, OIFITS_TOT_N)





    
#     #RAPPORTS D INTENSITE:    
    
#     distance_LM, intensity_LM, distance_N, intensity_N, wavel_LM, wavel_N = output_intensity_profile(PATH_OUTPUT_INTENSITY_FILE_LM, PATH_OUTPUT_INTENSITY_FILE_N, wavel_UD)
    
#     # LM BAND
#     fig1 = plt.figure()
#     W_LM,D_LM=np.meshgrid(np.append(wavel_LM,wavel_LM[-1]+np.diff(wavel_LM)[-1]),np.append([0],distance_LM))
#     plt.pcolor(W_LM,D_LM,np.array(intensity_LM),norm=colors.LogNorm(vmin=1E-4, vmax=intensity_LM.max()),shading='auto')
#     ax=plt.gca()
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*18, c='blue', linestyle='--')
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*28, c='blue',linestyle='--')
#     # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*42, c='red')
#     # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     plt.hlines(9,min(wavel_LM),max(wavel_LM),colors='red')
#     plt.hlines(r1rc*star_diam/2,min(wavel_LM),max(wavel_LM),colors='red')
#     plt.hlines(star_diam/2,min(wavel_LM),max(wavel_LM),colors='red',linestyles='--')
#     ax.set_ylim(0,max(distance_LM))
#     ax.set_ylabel('Angular distance from the star [mas]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µR]',fontsize = 11)
#     ax.minorticks_on()
#     rect=plt.Rectangle((3.05, 0), 0.25, 56,alpha=0.8, color='lightgrey')
#     ax.add_patch(rect)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     ax.set_ylim(1,max(distance_LM))
#     plt.colorbar()
#     plt.title("LM-Band 2D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'LM_intensity_ratio.jpeg')
#     # plt.close(fig1)
    
#     # N BAND
#     fig1=plt.figure()
#     W_N,D_N=np.meshgrid(np.append(wavel_N,wavel_N[-1]+np.diff(wavel_N)[-1]),np.append([0],distance_N))
#     plt.pcolor(W_N,D_N,np.array(intensity_N),norm=colors.LogNorm(vmin=1E-4, vmax=intensity_N.max()))
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*20, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*34, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*60, c='red')
#     # plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*95, c='red')
#     ax=plt.gca()
#     ax.set_ylim(1,max(distance_N))
#     ax.minorticks_on()
#     ax.set_ylabel('Angular distance from the star center [mas]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     plt.colorbar()
#     plt.title("N-Band 3D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'N_intensity_ratio.jpeg',bbox_inches='tight')
#     plt.close(fig1)


    
    
    
    
#     # #RAPPORTS D INTENSITE:
    
#     # WAVEL_LM, WAVEL_N, distance_wanted_LM, DUSTY_INTENSITY_LM, distance_wanted_N, DUSTY_INTENSITY_N = \
#     #     dusty_intensities(distance_LM,distance_N)    
    
#     # fig1=plt.figure()
#     # W,D=np.meshgrid(np.append(WAVEL_LM,WAVEL_LM[-1]+np.diff(WAVEL_LM)[-1]),np.append([0],distance_wanted_LM))
#     # plt.pcolor(W,D,np.array(DUSTY_INTENSITY_LM).T,norm=colors.LogNorm(vmin=1E-4, vmax=np.array(DUSTY_INTENSITY_LM).max()))
#     # ax=plt.gca()
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*22, c='red')
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*32, c='red')
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*42, c='red')
#     # # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     # ax.set_ylim(1,max(distance_wanted_LM))
#     # ax.set_ylabel('Angular Distance from the star [mas]',fontsize = 11)
#     # ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     # ax.minorticks_on()
#     # ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # plt.colorbar()
#     # plt.title("LM-Band DUSTY 2D Intensity Ratio")
#     # plt.savefig(PATH_OUTPUT_SPECTRA+'LM_intensity_ratio_DUSTY.jpeg')
#     # plt.close(fig1)
    
    
#     # fig1=plt.figure()
#     # W,D=np.meshgrid(np.append(WAVEL_N,WAVEL_N[-1]+np.diff(WAVEL_N)[-1]),np.append([0],distance_wanted_N))
#     # plt.pcolor(W,D,np.array(DUSTY_INTENSITY_N).T,norm=colors.LogNorm(vmin=1E-4, vmax=np.array(DUSTY_INTENSITY_N).max()))
#     # ax=plt.gca()
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*22, c='red')
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*32, c='red')
#     # # plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*42, c='red')
#     # # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     # ax.set_ylim(1,max(distance_wanted_N))
#     # ax.set_ylabel('Angular Distance from the star [mas]',fontsize = 11)
#     # ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     # ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # ax.minorticks_on()
#     # plt.colorbar()
#     # plt.title("N-Band DUSTY 2D Intensity Ratio")
#     # plt.savefig(PATH_OUTPUT_SPECTRA+'N_intensity_ratio_DUSTY.jpeg')
#     # plt.close(fig1)
    
    
    
    
#     # RESCALE 
#     new_intensity=[]
#     for i in range(len(wavel_LM)):
        
#         new_distance = distance_LM/wavel_LM[i]
        
#         # INTERP
        
#         intensity_lambda = interp1d(new_distance,intensity_LM[:,i],fill_value=0,bounds_error=False)
#         distance = np.linspace(np.amin(distance_LM/wavel_LM[0]),np.amax(distance_LM/wavel_LM[0]),100)
        
#         new_intensity.append(intensity_lambda(distance))
    
#     fig1=plt.figure()
#     W,D=np.meshgrid(np.append(wavel_LM,wavel_LM[-1]+np.diff(wavel_LM)[-1]),np.append([0],distance))
#     plt.pcolor(W,D,np.array(new_intensity).T,norm=colors.LogNorm(vmin=1E-4, vmax=np.array(new_intensity).max()))
#     ax=plt.gca()
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*22, c='red')
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*32, c='red')
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*42, c='red')
#     # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     ax.set_ylim(1,max(distance))
#     ax.set_ylabel('Angular Distance from the star center/ wavel [mas/µm]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     ax.minorticks_on()
#     plt.colorbar()
#     plt.title("LM-Band 2D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'LM_intensity_ratio_r.jpeg')
#     plt.close(fig1)
    
#     med_intensity = np.median(new_intensity,axis=0)
    
#     plt.figure()
#     plt.plot(med_intensity)
#     ax=plt.gca()
#     ax.set_yscale('log')
    
#     new_intensity_2=[]
#     for i in range(len(wavel_LM)):
    
#         new_intensity_2.append(new_intensity[i]-med_intensity)
    
    
#     fig1=plt.figure()
#     W,D=np.meshgrid(np.append(wavel_LM,wavel_LM[-1]+np.diff(wavel_LM)[-1]),np.append([0],distance))
#     plt.pcolor(W,D,np.array(new_intensity_2).T,vmin=-1E-2, vmax=1E-2)#,norm=colors.LogNorm(vmin=1E-8, vmax=np.array(new_intensity).max()))
#     ax=plt.gca()
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*22, c='red')
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*32, c='red')
#     plt.plot(wavel_LM,wavel_LM*1E-6/baseline_max_LM/(wavel_LM[0]*1E-6/baseline_max_LM)*42, c='red')
#     # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     ax.set_ylim(1,max(distance))
#     ax.set_ylabel('Angular Distance from the star center/ wavel [mas/µm]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     ax.minorticks_on()
#     plt.colorbar()
#     plt.title("LM-Band 2D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'LM_intensity_ratio_r_s.jpeg')
#     plt.close(fig1)
    
    
    
    
#     # N BAND
    
     
#     from cmcrameri import cm
#     mp   = cm.lajolla
#     rmap = mp.reversed()     
#     # RESCALE 
#     new_intensity=[]
#     for i in range(len(wavel_N)):
        
#         new_distance = distance_N/wavel_N[i]
        
#         # INTERP
        
#         intensity_lambda = interp1d(new_distance,intensity_N[:,i],fill_value=0,bounds_error=False)
#         distance = np.linspace(np.amin(distance_N/wavel_N[0]),np.amax(distance_N/wavel_N[0]),100)
        
#         new_intensity.append(intensity_lambda(distance))
    
#     fig1=plt.figure()
#     W,D=np.meshgrid(np.append(wavel_N,wavel_N[-1]+np.diff(wavel_N)[-1]),np.append([0],distance))
#     plt.pcolor(W,D,np.array(new_intensity).T,norm=colors.LogNorm(vmin=1E-4, vmax=np.array(new_intensity).max()), cmap = rmap )
#     ax=plt.gca()
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*22, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*32, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*42, c='red')
#     # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     ax.set_ylim(1,max(distance))
#     ax.set_ylabel('Angular Distance from the star/ wavel [mas/µm]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     ax.minorticks_on()
#     plt.colorbar()
#     plt.title("N-Band 2D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'N_intensity_ratio_r.jpeg')
#     plt.close(fig1)
    
#     med_intensity = np.median(new_intensity,axis=0)
    
#     plt.figure()
#     plt.plot(med_intensity)
#     ax=plt.gca()
#     ax.set_yscale('log')
    
#     new_intensity_2=[]
#     for i in range(len(wavel_N)):
    
#         new_intensity_2.append(new_intensity[i]-med_intensity)
    
    
#     fig1=plt.figure()
#     W,D=np.meshgrid(np.append(wavel_N,wavel_N[-1]+np.diff(wavel_N)[-1]),np.append([0],distance))
#     plt.pcolor(W,D,np.array(new_intensity_2).T,vmin=-1E-2, vmax=1E-2)#,norm=colors.LogNorm(vmin=1E-8, vmax=np.array(new_intensity).max()))
#     ax=plt.gca()
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*22, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*32, c='red')
#     plt.plot(wavel_N,wavel_N*1E-6/baseline_max_N/(wavel_N[0]*1E-6/baseline_max_N)*42, c='red')
#     # plt.hlines(68,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(35,min(wavel_LM),max(wavel_LM),colors='red')
#     # plt.hlines(20,min(wavel_LM),max(wavel_LM),colors='red')
#     ax.set_ylim(1,max(distance))
#     ax.set_ylabel('Angular Distance from the star/ wavel [mas/µm]',fontsize = 11)
#     ax.set_xlabel('Wavelength [µm]',fontsize = 11)
#     ax.minorticks_on()
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#     plt.colorbar()
#     plt.title("N-Band 2D Intensity Ratio")
#     plt.savefig(PATH_OUTPUT_SPECTRA+'N_intensity_ratio_r_s.jpeg')
#     plt.close(fig1)
    

# ###########################################################################################
# # CHECK PIERRE FORMULA:



# wavel_P = wavel_UD[0]    
# index_UD = nearest_index(np.unique(wavel_TOT),wavel_P*1E-6)
# cond= wavel_TOT == np.unique(wavel_TOT)[index_UD]
# q      = q_UD  

# q_UD_fitting      = np.array(q_TOT[cond])  

# if wavel_P <6:
#     delta_tetha = 2
#     header,data = READ_PRETTY_TABLE(PATH_OUTPUT_INTENSITY_FILE_LM,1)

# else:
#     delta_theta = 5
#     header,data = READ_PRETTY_TABLE(PATH_OUTPUT_INTENSITY_FILE_N,1)

# index_intensity = np.where(header==str(np.round(wavel_P,4)))[0]

# z = np.pi*mas_to_rad(delta_tetha)*q

# h_k_tmp = data[:,index_intensity]
# h_k = np.array(list(itertools.chain.from_iterable(h_k_tmp)))
# h_sum_spec = np.sum(np.array([h_k[k]*(2*k+1) for k in range(len(h_k))]))

# V_P = np.zeros(len(q))

# for i in range(len(q)):
#     a_k = 2*h_k/(z[i]*h_sum_spec)
#     V_P[i] = np.abs(np.sum([a_k[k]*((k+1)*j1(z[i]*(k+1)) - k*j1(z[i]*k)) for k in range(len(h_k))]))**2


    
# plt.figure()
# plt.plot(q,V_P)
# plt.plot(q,V_P)
# ax=plt.gca()
# ax.set_yscale('log')
# ax.set_ylim(1E-7,1.2)
# ax.set_xlim(0,5E7)


# # ESSAI POUR RECONSTRUIRE LES VISIBILITES MODELISES

# # EXTRAIRE LES RAPPORTS DE FLUX
# PATH_OUTPUT_FIT_LM = PATH_OUTPUT_tmp + NAME_FOLDER +'/FIT_RESULTS/fit_flux_ratio_LM.dat'
# PATH_OUTPUT_FIT_N  = PATH_OUTPUT_tmp + NAME_FOLDER +'/FIT_RESULTS/fit_flux_ratio_N.dat'

# header_flux_LM,data_flux_LM = READ_PRETTY_TABLE(PATH_OUTPUT_FIT_LM,1)
# header_flux_N,data_flux_N   = READ_PRETTY_TABLE(PATH_OUTPUT_FIT_N,1)

# wavel_flux_LM = data_flux_LM[:,0]
# wavel_flux_N = data_flux_N[:,0]

# flux_LM = data_flux_LM[:,1:]
# flux_N  = data_flux_N[:,1:]

# V_model = []
# wavel_model = []
# baseline_model = []

# for k in range(len(wavel_UD)):

#     index_UD = nearest_index(np.unique(wavel_TOT),wavel_UD[k]*1E-6)
#     cond= wavel_TOT == np.unique(wavel_TOT)[index_UD]
#     q_UD_fitting      = np.array(q_TOT[cond])  
#     baseline          = q_UD_fitting*wavel_UD[k]*1E-6
#     # print(index_UD)
#     if wavel_UD[k]<6:
#         index = k  
#         # print(index)
#         res_UD     = diam_outter_ring[0]
#         res_UD_bef = diam_inner_ring[0]
#         res_F      = flux_LM[index,:] 
#         nb_UD      = len(diam_outter_ring[0])
#         CF         = np.sum([BIG_V2[k][i]*res_F[i] for i in range(nb_UD)],axis=0)
#         Vis2       = (CF)**2
        
#     else :
#         index      = k-len(wavel_flux_LM)
#         res_UD     = diam_outter_ring[1]
#         res_UD_bef = diam_inner_ring[1]
#         res_F      = flux_N[index,:] 
#         nb_UD      = len(diam_outter_ring[1])
#         CF         = np.sum([BIG_V2[k][i]*res_F[i] for i in range(nb_UD)],axis=0)
#         Vis2       = (CF)**2
    
     
#     V_model.append(Vis2)
#     wavel_model.append([wavel_UD[k]*1E-6]*len(Vis2))    
#     baseline_model.append(q_UD_fitting*wavel_UD[k]*1E-6)    

# V_model = list(itertools.chain.from_iterable(V_model))
# wavel_model = list(itertools.chain.from_iterable(wavel_model))
# baseline_model = list(itertools.chain.from_iterable(baseline_model))

# V2_model_TOT = np.array([wavel_model,baseline_model,V_model], dtype="object").T   

# # CONDITION FOR MATISSE DATA
# from copy import deepcopy
# from p0_data_reading_new_data import OIFITS_READING
# V2_MATISSE, UV, UV_TP, TP_MATISSE, FLUX_MATISSE_LM,FLUX_MATISSE_N, BANDWIDTH_MATISSE, OIFITS_TOT_LM, OIFITS_TOT_N = OIFITS_READING(MATISSE_DIR_LM, MATISSE_DIR_N)

# V2_MATISSE_plot = deepcopy(V2_MATISSE)
# MATISSE_argsort = np.argsort(V2_MATISSE_plot[:,0])
# for i in range(np.shape((V2_MATISSE))[1]):
#     V2_MATISSE_plot[:,i] = V2_MATISSE_plot[MATISSE_argsort,i]


# V2_model_plot = deepcopy(V2_model_TOT)
# MODEL_argsort = np.argsort(V2_model_TOT[:,0])
# for i in range(np.shape((V2_model_TOT))[1]):
#     V2_model_plot[:,i] = V2_model_plot[MODEL_argsort,i]

# cond_LM = V2_model_plot[:,0]<5E-6
# cond_N  = V2_model_plot[:,0]>5E-6

# plt.figure()

# baseline_list= [10.6,30.9,90,130]
# colors_plot = ["orange",'ForestGreen','DeepSkyBlue',"red"]
# BAND= 'N'

# for k in range(len(baseline_list)):
#     if BAND == 'LM':
#         cond = cond_LM
#     else:
#         cond= cond_N
        
#     baseline =np.unique(V2_MATISSE_plot[cond,1])[nearest_index(np.unique(V2_MATISSE_plot[cond,1]),baseline_list[k])]
#     print(baseline)
#     cond_baseline_model = V2_model_plot[:,1]==baseline
#     cond_TOT_model = np.logical_and(cond_baseline_model,cond)    


#     cond_baseline_MATISSE = V2_MATISSE_plot[:,1]==baseline
#     cond_TOT_MATISSE = np.logical_and(cond_baseline_MATISSE,cond)    


#     if BAND=='LM':
    
#         plt.plot(V2_model_plot[cond_TOT_model,0]*1E6,V2_model_plot[cond_TOT_model,2],label='B=%i m'%(baseline),c=colors_plot[k])
#         plt.scatter(V2_MATISSE_plot[cond_TOT_MATISSE,0]*1E6, V2_MATISSE_plot[cond_TOT_MATISSE,2],c=colors_plot[k])
#         ax=plt.gca()
#         ax.set_yscale('log')
#         ax.set_ylim([1E-6,3E1])
#         ax.set_xlim([2.8,3.7])
#     else:

#         plt.plot(V2_model_plot[cond_TOT_model,0]*1E6,V2_model_plot[cond_TOT_model,2],label='B=%i m'%(baseline),c=colors_plot[k])
#         plt.errorbar(V2_MATISSE_plot[cond_TOT_MATISSE,0]*1E6, V2_MATISSE_plot[cond_TOT_MATISSE,2],fmt='o',yerr=V2_MATISSE_plot[cond_TOT_MATISSE,3],c=colors_plot[k])
#         ax=plt.gca()
#         ax.set_ylim([-0.1,1.1])
#         ax.set_xlim([9.0,12.0])
        
# ax.set_xlabel('Wavelength [µm]')
# ax.set_ylabel('V2')
# plt.legend()

# # plt.figure()
# # plt.plot(q_model_plot[np.argsort(q_model_plot)],V2_model_plot[np.argsort(q_model_plot)])
# # plt.scatter(q_UD_fitting, V2_MATISSE_plot,s=2)
    