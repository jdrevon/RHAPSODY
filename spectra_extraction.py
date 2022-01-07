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

# FLUX_LM = [OIFITS_TOT_LM[i]['FLUX'] for i in range(len(OIFITS_TOT_LM))]
# FLUX_N  = [OIFITS_TOT_N[i]['FLUX'] for i in range(len(OIFITS_TOT_N))]


# OIFITS_SORTING()
# V2_MATISSE, UV, UV_TP, TP_MATISSE, FLUX_MATISSE_LM,FLUX_MATISSE_N, BAND_WIDTH_MATISSE, OIFITS_TOT_LM, OIFITS_TOT_N = OIFITS_READING()


def spectra_matisse(OIFITS_TOT_LM,OIFITS_TOT_N, PATH_OUTPUT_FIT_RES, PATH_OUTPUT_SPECTRA):
    
    
    if ntelescope == 2:

        AT=['AT1','AT2']
        
    if ntelescope == 3:

        AT=['AT1','AT2','AT3']
        
    if ntelescope == 4:
    
        AT=['AT1','AT2','AT3','AT4']
    
    # What we want is to make statistics per AT on the flux to flag all the flux data
    # So we will stock all the data per AT make the median of all the flux in order to have N median spectras for the N AT's 
    
    
    # We start looping over the AT's
    
    FLUX_LM      = np.zeros((len(OIFITS_TOT_LM),ntelescope,len(OIFITS_TOT_LM[0]['WAVEL'])))
    FLUX_LM_ERR      = np.zeros((len(OIFITS_TOT_LM),ntelescope,len(OIFITS_TOT_LM[0]['WAVEL'])))
    
    FLUX_N       = np.zeros((len(OIFITS_TOT_N),ntelescope,len(OIFITS_TOT_N[0]['WAVEL'])))
    FLUX_N_ERR       = np.zeros((len(OIFITS_TOT_N),ntelescope,len(OIFITS_TOT_N[0]['WAVEL'])))
    
    wavel_LM      = np.zeros((len(OIFITS_TOT_LM),ntelescope,len(OIFITS_TOT_LM[0]['WAVEL'])))
    wavel_N       = np.zeros((len(OIFITS_TOT_N),ntelescope,len(OIFITS_TOT_N[0]['WAVEL'])))
    
    j=0

    for item in AT:
         
        for i in range(len(OIFITS_TOT_LM)):
            
            # I select the flux data from the i-th AT
            AT_index = np.where(OIFITS_TOT_LM[i]['FLUX']['AT_NUMBER']==item)[0][0]
            
            # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
            # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
            
            FLAG_LM = OIFITS_TOT_LM[i]['FLUX']['FLAG'][AT_index]
            
            FLUX_LM[i][j][:]  = ma.array(OIFITS_TOT_LM[i]['FLUX']['FLUX'][AT_index], mask= FLAG_LM)
            FLUX_LM_ERR[i][j][:]  = ma.array(OIFITS_TOT_LM[i]['FLUX']['FLUX_ERR'][AT_index], mask= FLAG_LM)
            
            wavel_LM[i][j][:] = ma.array(OIFITS_TOT_LM[i]['FLUX']['WAVEL'][AT_index], mask= FLAG_LM)
                
        for i in range(len(OIFITS_TOT_N)):
            
            # Same as above for N-band

            # I select the flux data from the i-th AT
            AT_index = np.where(OIFITS_TOT_N[i]['FLUX']['AT_NUMBER']==item)[0][0]
            
            # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
            # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
            FLAG_N = OIFITS_TOT_N[i]['FLUX']['FLAG'][AT_index]
                        
            FLUX_N[i][j][:]  = ma.array(OIFITS_TOT_N[i]['FLUX']['FLUX'][AT_index], mask= FLAG_N)
            FLUX_N_ERR[i][j][:]  = ma.array(OIFITS_TOT_N[i]['FLUX']['FLUX_ERR'][AT_index], mask= FLAG_N)
            
            wavel_N[i][j][:] = ma.array(OIFITS_TOT_N[i]['FLUX']['WAVEL'][AT_index], mask= FLAG_N)
        

        # FOR LM BAND:
         
        x_i_LM = ma.array(FLUX_LM)    
         
        sigma_sup_LM = ma.array(FLUX_LM_ERR)
        W_LM = 1/sigma_sup_LM**2 
        X_w_LM = np.sum(W_LM*x_i_LM,axis=0)/np.sum(W_LM,axis=0)
        V_X_w_LM = np.sum(W_LM**2,axis=0)/(np.sum(W_LM,axis=0)**2-np.sum(W_LM**2,axis=0))*np.sum(W_LM*(x_i_LM-X_w_LM)**2,axis=0)/np.sum(W_LM,axis=0) 
        sigma_LM = np.sqrt(V_X_w_LM)
        
        # FOR N BAND:
        
        x_i_N = ma.array(FLUX_N)    
        
        sigma_sup_N = ma.array(FLUX_N_ERR)
        W_N = 1/sigma_sup_N**2 
        X_w_N = np.sum(W_N*x_i_N,axis=0)/np.sum(W_N,axis=0)
        V_X_w_N = np.sum(W_N**2,axis=0)/(np.sum(W_N,axis=0)**2-np.sum(W_N**2,axis=0))*np.sum(W_N*(x_i_N-X_w_N)**2,axis=0)/np.sum(W_N,axis=0) 
        sigma_N = np.sqrt(V_X_w_N)


        # FLUX_N_median  = ma.median(ma.array(FLUX_N),axis=0)
        # FLUX_LM_median = ma.median(ma.array(FLUX_LM),axis=0)
    
        # FLUX_N_std  = ma.median(ma.abs(ma.array(FLUX_N)-FLUX_N_median),axis=0)*1.48
        # FLUX_LM_std = ma.median(ma.abs(ma.array(FLUX_LM)-FLUX_LM_median),axis=0)*1.48
        

        WAVEL         = [*(OIFITS_TOT_LM[0]['WAVEL']*1E6),*(OIFITS_TOT_N[0]['WAVEL']*1E6)]
        FLUX_DATA     = [*X_w_LM[j],*X_w_N[j]]
        FLUX_DATA_err = [*sigma_LM[j],*sigma_N[j]]


        plt.figure()
        plt.errorbar(WAVEL,FLUX_DATA,yerr=FLUX_DATA_err, fmt='o', ms=1)
        ax=plt.gca()
        ax.set_xlabel('Wavelengths [µm]')
        ax.set_ylabel('Flux [Jy]')
        ax.set_xscale('log')
        ax.set_ylim(bottom=-20)
        plt.savefig(PATH_OUTPUT_SPECTRA+'%s_spectra.jpg'%item,bbox_inches='tight', dpi=400)
        plt.close()


        data_flux = PrettyTable()
        data_flux.field_names = ["Wvl [µm]"]+['Flux [Jy]']+['ERR Flux [Jy]']
        for k in range(len(WAVEL)):
            data_flux.add_row([WAVEL[k]]+ [FLUX_DATA[k]] + [FLUX_DATA_err[k]])
        with open(PATH_OUTPUT_FIT_RES+'%s_FLUX.dat'%item, 'w') as f: f.write(str(data_flux))


        j=j+1


    FLUX_N_TOT  = ma.mean(ma.array(X_w_N),axis=0)
    FLUX_LM_TOT = ma.mean(ma.array(X_w_LM),axis=0)

    FLUX_N_std_TOT  = ma.mean(ma.array(sigma_LM), axis=0)
    FLUX_LM_std_TOT = ma.mean(ma.array(sigma_N), axis=0)

    WAVEL         = [*(OIFITS_TOT_LM[0]['WAVEL']*1E6),*(OIFITS_TOT_N[0]['WAVEL']*1E6)]
    FLUX_DATA     = [*FLUX_LM_TOT,*FLUX_N_TOT]
    FLUX_DATA_err = [*FLUX_LM_std_TOT,*FLUX_N_std_TOT]

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

