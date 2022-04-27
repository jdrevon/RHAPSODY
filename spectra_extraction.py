#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 09:09:48 2021

@author: jdrevon
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from rhapsody_init import ntelescope
from prettytable import PrettyTable


def spectra_data(OIFITS_TOT, PATH_OUTPUT_FIT_RES, PATH_OUTPUT_SPECTRA):
    
    
    if ntelescope == 2:

        AT=['AT1','AT2']
        
    if ntelescope == 3:

        AT=['AT1','AT2','AT3']
        
    if ntelescope == 4:
    
        AT=['AT1','AT2','AT3','AT4']
    
    # What we want is to make statistics per AT on the flux to flag all the flux data
    # So we will stock all the data per AT make the median of all the flux in order to have N median spectras for the N AT's 
    
    
    # We start looping over the AT's
    
    FLUX           = np.zeros((len(OIFITS_TOT),ntelescope,len(OIFITS_TOT[0]['WAVEL'])))
    FLUX_ERR       = np.zeros((len(OIFITS_TOT),ntelescope,len(OIFITS_TOT[0]['WAVEL'])))
    wavel          = np.zeros((len(OIFITS_TOT),ntelescope,len(OIFITS_TOT[0]['WAVEL'])))
    

    for j in range(len(AT)):
         
        for i in range(len(OIFITS_TOT)):
            
            # I select the flux data from the i-th AT
            AT_index = np.where(OIFITS_TOT[i]['FLUX']['AT_NUMBER']==AT[j])[0][0]
            
            # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX[i][j][:]
            # So in the first loop AT1 values will be stocked in the first slot of the FLUX array.
            
            FLAG = OIFITS_TOT[i]['FLUX']['FLAG'][AT_index]
            
            FLUX[i][j]      = ma.array(OIFITS_TOT[i]['FLUX']['FLUX'][AT_index], mask= FLAG)
            FLUX_ERR[i][j]  = ma.array(OIFITS_TOT[i]['FLUX']['FLUX_ERR'][AT_index], mask= FLAG)
            wavel[i][j]     = ma.array(OIFITS_TOT[i]['FLUX']['WAVEL'][AT_index], mask= FLAG)
        
        # Each AT spectra are the weighted mean of all the indivudal spectra with this telescope
         
        x_i = ma.array(FLUX)    
         
        sigma_sup = ma.array(FLUX_ERR)
        W = 1/sigma_sup**2 
        X_w = np.sum(W*x_i,axis=0)/np.sum(W,axis=0)
        V_X_w = np.sum(W**2,axis=0)/(np.sum(W,axis=0)**2-np.sum(W**2,axis=0))*np.sum(W*(x_i-X_w)**2,axis=0)/np.sum(W,axis=0) 
        sigma = np.sqrt(V_X_w)
                

        WAVEL         = OIFITS_TOT[0]['WAVEL']*1E6
        FLUX_DATA     = X_w[j]
        FLUX_DATA_err = sigma[j]


        plt.figure()
        plt.errorbar(WAVEL,FLUX_DATA,yerr=FLUX_DATA_err, fmt='o', ms=1)
        ax=plt.gca()
        ax.set_xlabel('Wavelengths [µm]')
        ax.set_ylabel('Flux [Jy]')
        ax.set_xscale('log')
        ax.set_ylim(bottom=-20)
        plt.savefig(PATH_OUTPUT_SPECTRA+'%s_spectra.jpg'%AT[j],bbox_inches='tight', dpi=400)
        plt.close()


        data_flux = PrettyTable()
        data_flux.field_names = ["Wvl [µm]"]+['Flux [Jy]']+['ERR Flux [Jy]']
        for k in range(len(WAVEL)):
            data_flux.add_row([WAVEL[k]]+ [FLUX_DATA[k]] + [FLUX_DATA_err[k]])
        with open(PATH_OUTPUT_FIT_RES+'%s_FLUX.dat'%AT[j], 'w') as f: f.write(str(data_flux))



    FLUX_TOT      = ma.mean(ma.array(X_w),axis=0)
    FLUX_std_TOT  = ma.mean(ma.array(sigma), axis=0)

    # The total spectra is the mean of all the AT's spectra computed before

    FLUX_DATA     = FLUX_TOT
    FLUX_DATA_err = FLUX_std_TOT

    data_flux = PrettyTable()
    data_flux.field_names = ["Wvl [µm]"]+['Flux [Jy]']+['ERR Flux [Jy]']
    for k in range(len(WAVEL)):
        data_flux.add_row([WAVEL[k]]+ [FLUX_DATA[k]] + [FLUX_DATA_err[k]])
    with open(PATH_OUTPUT_FIT_RES+'MEAN_FLUX.dat', 'w') as f: f.write(str(data_flux))

    plt.figure()
    plt.errorbar(WAVEL,FLUX_DATA,yerr=FLUX_DATA_err, fmt='o', ms=1)
    ax=plt.gca()
    ax.set_xlabel('Wavelengths [µm]')
    ax.set_ylabel('Flux [Jy]')
    ax.set_xscale('log')
    ax.set_ylim(bottom=-20)
    plt.savefig(PATH_OUTPUT_SPECTRA+'AT_mean_spectra.jpg',bbox_inches='tight', dpi=400)
    plt.close()


    return WAVEL, FLUX_DATA, FLUX_DATA_err

