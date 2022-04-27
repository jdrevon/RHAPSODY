#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 18:15:26 2021

@author: jdrevon
"""

import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from rhapsody_init import DATA_band_name

def plot_intensity_model_image(x_model, image_model, wavel_UD, R_image, **kwargs):

    'help(function_name) for further explanation\n \
        kwargs\n \
        /PLOT\ Show the plot or not [default = False]\n \
        /SAVE_OUTPUT\  Full path for figures saving. [default = None]\n \
        '    

    PLOT        = kwargs.get('PLOT', False)    
    SAVE_OUTPUT = kwargs.get('SAVE_OUTPUT', None)    
        
       
    fig1, axs = plt.subplots()

    im=axs.imshow(abs(image_model),extent=[min(x_model),max(x_model),min(x_model),max(x_model)], cmap = 'hot',norm=colors.PowerNorm(gamma=0.5))#, norm=colors.PowerNorm(gamma=1.5))
    cb = fig1.colorbar(im) 
    cb.set_label(label='Intensity Ratio [I_tot/I_max]', size=16, labelpad=10)
    cb.ax.tick_params(labelsize= 16)

    axs.set_xlabel(r'$\alpha$ [mas]',fontsize = 16)
    axs.set_ylabel(r'$\delta$ [mas]',fontsize = 16)
    axs.set_title(r'Image from the model at $\lambda$ = %.2f µm'%wavel_UD,fontsize = 16)

    axs.minorticks_on()
    axs.tick_params(axis='x', labelsize=13)
    axs.tick_params(axis='y', labelsize=13)
    axs.set_xlim([-R_image,R_image])
    axs.set_ylim([-R_image,R_image])

    if SAVE_OUTPUT != None:
        
        fig1.savefig(SAVE_OUTPUT +'image' + str(np.round(wavel_UD,3)) + '.jpg', bbox_inches = 'tight')
   
    if PLOT == False:
        plt.close(fig1) 
   
    return




def plot_intensity_model_profile(x_model, y_model, wavel_UD, **kwargs):

    'help(function_name) for further explanation\n \
        kwargs\n \
        /PLOT\ Show the plot or not [default = False]\n \
        /SAVE_OUTPUT\  Full path for figures saving. [default = None]\n \
        /logx\  Logarithmic scale on the x-axis. [default = False]\n \
        /logy\  Logarithmic scale on the y-axis. [default = False]\n \
        /x_rec_err\  Array of the error on the observations for the x-axis. [default = None]\n \
        /y_rec_err\  Array of the error on the observations for the y-axis. [default = None]\n \
        /xlim_min\  Minimum x_value of the plot window [default = min([x_rec,x_model])]\n \
        /xlim_max\  Maximum x_value of the plot window [default = max([x_rec,x_model])]\n \
        /ylim_min\  Minimum y_value of the plot window [default = min([y_rec,y_model])]\n \
        /ylim_max\  Maximum y_value of the plot window [default = max([y_rec,y_model])]\n \
        '    

    PLOT        = kwargs.get('PLOT', False)    
    SAVE_OUTPUT = kwargs.get('SAVE_OUTPUT', None)    
    logx        = kwargs.get('logx', False)
    logy        = kwargs.get('logy', False)
    xlim_min    = kwargs.get('xlim_min', -0.05)
    xlim_max    = kwargs.get('xlim_max', None)
    ylim_min    = kwargs.get('ylim_min', -0.05)
    ylim_max    = kwargs.get('ylim_max', 1.05)
        
        
    fig2 = plt.figure()
    axs2 = plt.gca()
    axs2.plot(x_model,y_model, label='Reconstructive Profile')
    axs2.set_xlabel(r' Angular distance [mas]',fontsize = 16)
    axs2.set_ylabel(r'$I_{im}/I_{max}$',fontsize = 16)
    axs2.set_title(r'Radial Profile at $\lambda$ = %.2f µm'%wavel_UD,fontsize = 16)
    axs2.minorticks_on()
    axs2.tick_params(axis='x', labelsize=13)
    axs2.tick_params(axis='y', labelsize=13)
    axs2.set_xlim([xlim_min,xlim_max])
    axs2.set_ylim([ylim_min,ylim_max])
    axs2.legend(loc='upper right')
    
    if logx==True:
        axs2.set_xscale('log')

    if logy==True:
        axs2.set_yscale('log')

    
    if SAVE_OUTPUT != None:

        fig2.savefig(SAVE_OUTPUT +'radial_' + str(np.round(wavel_UD,3)) + '.jpg', bbox_inches = 'tight')
   
    if PLOT == False:
        plt.close(fig2) 
   
    return



