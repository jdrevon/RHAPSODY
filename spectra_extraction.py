# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 09:09:48 2021

@author: jdrevon
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from astropy.table import Table
import astropy.io.ascii as ascii_f

# FLUX_LM = [OIFITS_TOT_LM[i]['FLUX'] for i in range(len(OIFITS_TOT_LM))]
# FLUX_N  = [OIFITS_TOT_N[i]['FLUX'] for i in range(len(OIFITS_TOT_N))]


# OIFITS_SORTING()
# V2_MATISSE, UV, UV_TP, TP_MATISSE, FLUX_MATISSE_LM,FLUX_MATISSE_N, BAND_WIDTH_MATISSE, OIFITS_TOT_LM, OIFITS_TOT_N = OIFITS_READING()


def spectra_matisse(OIFITS_TOT_LM,OIFITS_TOT_N):
    
    FLUX_LM          = []
    FLUX_LM_err      = []
    
    FLUX_N           = []
    FLUX_N_err       = []
    
    wavel_LM         = []
    wavel_N          = []
    
    for i in range(len(OIFITS_TOT_LM)):
        
        cond = OIFITS_TOT_LM[i]['FLUX']['FLAG']==False
        for k in range(len(cond)):
            FLUX_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['FLUX'][k].astype('float'),mask=~cond[k]))
            FLUX_LM_err.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['FLUX_ERR'][k].astype('float'),mask=~cond[k]))
            wavel_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['WAVEL'][k].astype('float'),mask=~cond[k]))
    
    
    for i in range(len(OIFITS_TOT_N)):
        
        cond = OIFITS_TOT_N[i]['FLUX']['FLAG']==False
        for k in range(len(cond)):
            FLUX_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['FLUX'][k].astype('float'),mask=~cond[k]))
            FLUX_N_err.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['FLUX_ERR'][k].astype('float'),mask=~cond[k]))
            wavel_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['WAVEL'][k].astype('float'),mask=~cond[k]))
    
    # FOR LM BAND:
    
    x_i_LM = ma.array(FLUX_LM)    
    # X_LM   = np.sum(x_i_LM,axis=0)/len(x_i_LM)
    # V_X_LM = np.sum((x_i_LM-X_LM)**2,axis=0)/len(x_i_LM)/(len(x_i_LM)-1)
    
    sigma_sup_LM = ma.array(FLUX_LM_err)
    W_LM = 1/sigma_sup_LM**2 #np.array([1]*len(x_i))#1/sigma_sup**2 #np.array([1]*len(x_i)) # #
    X_w_LM = np.sum(W_LM*x_i_LM,axis=0)/np.sum(W_LM,axis=0)
    V_X_w_LM = np.sum(W_LM**2,axis=0)/(np.sum(W_LM,axis=0)**2-np.sum(W_LM**2,axis=0))*np.sum(W_LM*(x_i_LM-X_w_LM)**2,axis=0)/np.sum(W_LM,axis=0) 
    sigma_LM = np.sqrt(V_X_w_LM)
    
    # FOR N BAND:
    
    x_i_N = ma.array(FLUX_N)    
    # X_N   = np.sum(x_i_N,axis=0)/len(x_i_N)
    # V_X_N = np.sum((x_i_N-X_N)**2,axis=0)/len(x_i_N)/(len(x_i_N)-1)
    
    sigma_sup_N = ma.array(FLUX_N_err)
    W_N = 1/sigma_sup_N**2 #np.array([1]*len(x_i))#1/sigma_sup**2 #np.array([1]*len(x_i)) # #
    X_w_N = np.sum(W_N*x_i_N,axis=0)/np.sum(W_N,axis=0)
    V_X_w_N = np.sum(W_N**2,axis=0)/(np.sum(W_N,axis=0)**2-np.sum(W_N**2,axis=0))*np.sum(W_N*(x_i_N-X_w_N)**2,axis=0)/np.sum(W_N,axis=0) 
    sigma_N = np.sqrt(V_X_w_N)
    

    # SAVING DATA:
    
    WAVEL         = [*(OIFITS_TOT_LM[0]['WAVEL']*1E6),*(OIFITS_TOT_N[0]['WAVEL']*1E6)]
    FLUX_DATA     = [*X_w_LM,*X_w_N]
    FLUX_DATA_err = [*sigma_LM,*sigma_N]
    
    data = Table([WAVEL, FLUX_DATA, FLUX_DATA_err], names=['WAVEL [µm]','FLUX [Jy]', 'FLUX_err [Jy]'])
    ascii_f.write(data, 'FLUX_MATISSE_R_Scl.dat', overwrite=True)

    return WAVEL, FLUX_DATA, FLUX_DATA_err

