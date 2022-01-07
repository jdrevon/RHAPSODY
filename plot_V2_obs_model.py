#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:40:12 2021

@author: jdrevon
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_V2_obs_model(x_obs,y_obs, x_model, y_model, wavel_obs, **kwargs):

    'help(function_name) for further explanation\n \
        kwargs\n \
        /PLOT\ Show the plot or not [default = False]\n \
        /SAVE_OUTPUT\  Full path for figures saving. [default = None]\n \
        /logx\  Logarithmic scale on the x-axis. [default = False]\n \
        /logy\  Logarithmic scale on the y-axis. [default = False]\n \
        /x_obs_err\  Array of the error on the observations for the x-axis. [default = None]\n \
        /y_obs_err\  Array of the error on the observations for the y-axis. [default = None]\n \
        /xlim_min\  Minimum x_value of the plot window [default = min([x_obs,x_model])]\n \
        /xlim_max\  Maximum x_value of the plot window [default = max([x_obs,x_model])]\n \
        /ylim_min\  Minimum y_value of the plot window [default = min([y_obs,y_model])]\n \
        /ylim_max\  Maximum y_value of the plot window [default = max([y_obs,y_model])]\n \
        '    

    PLOT        = kwargs.get('PLOT', False)    
    SAVE_OUTPUT = kwargs.get('SAVE_OUTPUT', None)    
    logx        = kwargs.get('logx', False)
    logy        = kwargs.get('logy', False)
    x_obs_err   = kwargs.get('x_obs_err', None)
    y_obs_err   = kwargs.get('y_obs_err', None)
    xlim_min    = kwargs.get('xlim_min', min(np.concatenate([x_obs,x_model])))
    xlim_max    = kwargs.get('xlim_max', max(np.concatenate([x_obs,x_model])))
    ylim_min    = kwargs.get('ylim_min', min(np.concatenate([y_obs,y_model])))
    ylim_max    = kwargs.get('ylim_max', max(np.concatenate([y_obs,y_model])))
    
    
    # print(PLOT,SAVE_OUTPUT,logx,logy,xlim_min,xlim_max,ylim_min, ylim_max)
    
    
    fig1 = plt.figure()
    axs1 = plt.gca()
    
    axs1.plot(x_model,y_model, label='Model')
    axs1.errorbar(x_obs,y_obs, xerr = x_obs_err, yerr= y_obs_err, fmt='o', ms=2, ecolor='lightgray', c='orange', label="Observations")#, norm=colors.PowerNorm(gamma=1.5))
    axs1.set_xlabel(r' $B/\lambda$ [rad$^{-1}$]',fontsize = 16)
    axs1.set_ylabel(r'V2',fontsize = 16)

    if logx==True:
        axs1.set_xscale('log')

    if logy==True:
        axs1.set_yscale('log')

    axs1.set_xlim([xlim_min,xlim_max])
    axs1.set_ylim([ylim_min,ylim_max])
    axs1.set_title(r'V2 at $\lambda$ = %.2f µm'%(np.round(wavel_obs,3)),fontsize = 16)
    axs1.minorticks_on()
    axs1.tick_params(axis='x', labelsize=13)
    axs1.tick_params(axis='y', labelsize=13)
    axs1.legend(loc='upper right')
    
    if SAVE_OUTPUT != None:
        
        if wavel_obs<6:
            if np.logical_or(logx,logy):            
                fig1.savefig(SAVE_OUTPUT +'V2_log_LM_' + str(np.round(wavel_obs,3)) + '.jpg', bbox_inches = 'tight')
            else:
                fig1.savefig(SAVE_OUTPUT +'V2_linear_LM_' + str(np.round(wavel_obs,3)) + '.jpg', bbox_inches = 'tight')
                
        else : 
            if np.logical_or(logx,logy):            
                fig1.savefig(SAVE_OUTPUT +'V2_log_N_' + str(np.round(wavel_obs,3)) + '.jpg', bbox_inches = 'tight')
            else:
                fig1.savefig(SAVE_OUTPUT +'V2_linear_N_' + str(np.round(wavel_obs,3)) + '.jpg', bbox_inches = 'tight')
            
    
    if PLOT==False:
        plt.close(fig1)
        
    return