#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 01:24:13 2021

@author: jdrevon
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_histogram(V2_obs,V2_obs_err,V2_model, wavel_UD, **kwargs):
    
    PLOT        = kwargs.get('PLOT', False)    
    SAVE_OUTPUT = kwargs.get('SAVE_OUTPUT', None)    

    
    fig2 = plt.figure()
    V2_hist=(V2_model-V2_obs)/V2_obs_err  
    median_hist = np.median(V2_hist)
    MAD = np.median(np.abs(V2_hist-median_hist))               
    # plt.hist(V2_hist[cond2]/V2_obs_err[cond2],bins=(10**np.linspace(-3, 2, 100)),label='weighted residuals')
    bins_nb= int(len(V2_hist)/15)
    plt.hist(V2_hist,bins=bins_nb,label='Pearson residuals')
    plt.vlines(median_hist,ymin = 0, ymax = len(V2_obs_err)/2, label='Med = %.2e'%median_hist, linestyles='dashed', alpha = 0.8)
    plt.vlines(median_hist-4.5*MAD,ymin = 0, ymax = len(V2_obs_err)/2,label='-4.5*MAD', linestyles='dashed', colors='red',alpha=0.4)
    plt.vlines(median_hist+4.5*MAD,ymin = 0, ymax = len(V2_obs_err)/2,label='+4.5*MAD',linestyles='dashed',colors='red',alpha=0.4)
    # plt.xscale('log')
    ax=plt.gca()
    ax.set_xlabel('Pearson residuals')
    ax.set_xlim(median_hist-7*MAD,median_hist+7*MAD)
    plt.title(r'$\lambda$ = %0.3f µm'%(np.round(wavel_UD,3)))
    plt.legend()

    if PLOT == False:
        plt.close(fig2)

    if SAVE_OUTPUT != None:
    
        fig2.savefig(SAVE_OUTPUT +'hist' + str(np.round(wavel_UD,3)) + '.jpg', bbox_inches = 'tight')
            
