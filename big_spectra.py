# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:50:08 2022

@author: jdrevon
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib
import matplotlib.colors as colors
from cmcrameri import cm
import numpy as np
from black_body import BB_v_Jy
from rhapsody_init import BB_norm, BB_temperature, stellar_radii, distance_target
from A1_mas_to_rad import au_to_R_sun

def big_spectra_with_flux(W,D,intensity, FLUX_WAVEL, FLUX_DATA, PATH_OUTPUT_SPECTRA):
    
    widths  = [1]
    heights = [1,5]

    mp   = cm.lajolla
    rmap = mp.reversed()     

    fig=plt.figure(figsize=(4, 12))
    gs1 = GridSpec(2, 1, width_ratios=widths, height_ratios=heights)
        

    ax3 = plt.subplot(gs1[0,0])
    ax1 = plt.subplot(gs1[1,0])    

    # First part of the plot deal with intensity profile
    
    plot1=ax1.pcolor(W,D,intensity,norm=colors.LogNorm(vmin=1E-4, vmax=1),shading='auto', cmap=rmap)
            
    ax1.set_yscale('log')
    ax1.set_ylim(3,np.amax(D))
    # ax1.set_ylim(3,100)
    ax1.yaxis.set_ticks_position('both')




    ax1.set_ylabel('Angular distance from the star center [mas]',fontsize = 13)
    ax1.set_yticks([5,10,12,15,18,25,30,40,50,75,100])    
    ax1.tick_params(axis='both', labelsize=13)
    ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.05, 0.7])
    fig.colorbar(plot1, cax=cbar_ax).set_label(label=r'$I(\rho)/I(0)$',size=13)
    fig.colorbar(plot1, cax=cbar_ax).ax.tick_params(labelsize=13) 


    wavel_MATISSE_flux = np.array(FLUX_WAVEL)    
    

    if BB_norm == True:
    
        wavel_BB, blackbody   = BB_v_Jy(wavel_MATISSE_flux*1E-6, BB_temperature, au_to_R_sun(stellar_radii), distance_target)
        data2                 = FLUX_DATA/blackbody.astype('float')
        MATISSE_flux       = data2/np.mean(data2)

    else:

        data2              = FLUX_DATA            
        MATISSE_flux       = data2/np.mean(data2)

    
    ax3.tick_params(axis='both', labelsize=15)
    ax3.set_xlim(np.amin(W),np.amax(W))
    ax3.set_ylim(0.5,1.5)
    ax3.set_ylabel('Normalized Flux',fontsize =15)

    ax3.plot(np.sort(wavel_MATISSE_flux),[x for _, x in sorted(zip(wavel_MATISSE_flux, MATISSE_flux))], label = 'MATISSE')
    ax3.minorticks_on()
    ax3.set

    handles, labels = [(a + b + c) for a, b, c in zip(ax1.get_legend_handles_labels(), ax1.get_legend_handles_labels(), ax3.get_legend_handles_labels())]
    fig1=ax1.figure
    fig1.text(0.5,0.04, 'Wavelength [µm]',fontsize = 15, ha="center", va="center")
    
    # plt.setp(ax1.get_yticklabels(), visible=False)

    plt.savefig(PATH_OUTPUT_SPECTRA+'big_spectra.jpg',bbox_inches='tight', dpi=400)
    plt.close(fig)
    
    return


def big_spectra_without_flux(W,D,intensity, PATH_OUTPUT_SPECTRA):
    

    mp   = cm.lajolla
    rmap = mp.reversed()     

    fig=plt.figure(figsize=(4, 12))

    ax1 = plt.gca()    

    # First part of the plot deal with intensity profile
    
    plot1=ax1.pcolor(W,D,intensity,norm=colors.LogNorm(vmin=1E-4, vmax=1),shading='auto', cmap=rmap)
    # plt.plot(W[0],W[0]/(W[0][0])*22, c='white', linestyle='--')
    plt.plot(W[0],W[0]/(W[0][0])*32, c='white',linestyle='--')
    plt.plot(W[0],W[0]/(W[0][0])*42, c='white',linestyle='--')
    plt.plot(W,W*1E-6/140/(W[0]*1E-6/140)*32, c='white',linestyle='--')
    plt.plot(W,W*1E-6/140/(W[0]*1E-6/140)*42, c='white',)

    ax1.axhline(y=12, color='chartreuse', linestyle='-')
    ax1.axhline(y=18, color='chartreuse', linestyle='-')

    ax1.axhline(y=24, color='deepskyblue', linestyle='-')
        
    
    ax1.set_yscale('log')
    ax1.set_ylim(3,np.amax(D))
    # ax1.set_ylim(3,100)
    ax1.yaxis.set_ticks_position('both')

    ax1.set_ylabel('Angular distance from the star center [mas]',fontsize = 13)
    ax1.set_yticks([5,10,12,15,18,25,30,40,50,75,100])    
    ax1.tick_params(axis='both', labelsize=13)
    ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())


    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.05, 0.7])
    fig.colorbar(plot1, cax=cbar_ax).set_label(label=r'$I(\rho)/I(0)$',size=13)
    fig.colorbar(plot1, cax=cbar_ax).ax.tick_params(labelsize=13) 

    fig.text(0.5,0.04, 'Wavelength [µm]',fontsize = 15, ha="center", va="center")

    # plt.setp(ax1.get_yticklabels(), visible=False)

    plt.savefig(PATH_OUTPUT_SPECTRA+'big_spectra.jpg',bbox_inches='tight', dpi=400)
    plt.close(fig)
    
    return
