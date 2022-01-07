#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:34:08 2020

@author: jdrevon
"""
import numpy as np
from rhapsody_init import ntelescope


def FLUX_FLAG(OIFITS_TOT_LM, OIFITS_TOT_N):
    
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
    FLUX_N       = np.zeros((len(OIFITS_TOT_N),ntelescope,len(OIFITS_TOT_N[0]['WAVEL'])))
    
    wavel_LM      = np.zeros((len(OIFITS_TOT_LM),ntelescope,len(OIFITS_TOT_LM[0]['WAVEL'])))
    wavel_N       = np.zeros((len(OIFITS_TOT_N),ntelescope,len(OIFITS_TOT_N[0]['WAVEL'])))
    
    j=0

    for item in AT:
         
        for i in range(len(OIFITS_TOT_LM)):
            
            # I select the flux data from the i-th AT
            AT_index = np.where(OIFITS_TOT_LM[i]['FLUX']['AT_NUMBER']==item)[0][0]
            
            # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
            # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
                        
            FLUX_LM[i][j][:]  = OIFITS_TOT_LM[i]['FLUX']['FLUX'][AT_index]
            
            wavel_LM[i][j][:] = OIFITS_TOT_LM[i]['FLUX']['WAVEL'][AT_index]
                
        for i in range(len(OIFITS_TOT_N)):
            
            # Same as above for N-band

            # I select the flux data from the i-th AT
            AT_index = np.where(OIFITS_TOT_N[i]['FLUX']['AT_NUMBER']==item)[0][0]
            
            # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
            # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
                        
            FLUX_N[i][j][:]  = OIFITS_TOT_N[i]['FLUX']['FLUX'][AT_index]
            
            wavel_N[i][j][:] = OIFITS_TOT_N[i]['FLUX']['WAVEL'][AT_index]           
        
        FLUX_N_median  = np.median(np.array(FLUX_N),axis=0)
        FLUX_LM_median = np.median(np.array(FLUX_LM),axis=0)
    
        FLUX_N_std  = np.median(np.abs(np.array(FLUX_N)-FLUX_N_median),axis=0)*1.48
        FLUX_LM_std = np.median(np.abs(np.array(FLUX_LM)-FLUX_LM_median),axis=0)*1.48


        print('SORTING LM BAND FOR '+ item)
        
        for file in OIFITS_TOT_LM:

            AT_index = np.where(OIFITS_TOT_LM[i]['FLUX']['AT_NUMBER']==item)[0][0]

            cond = np.logical_or(file['FLUX']['FLUX'][AT_index]>=FLUX_LM_median[j]+2*FLUX_LM_std[j],\
                                  file['FLUX']['FLUX'][AT_index]<=FLUX_LM_median[j]-2*FLUX_LM_std[j])


            file['FLUX']['FLAG'][AT_index] = cond | file['FLUX']['FLAG'][AT_index]
                

        print('END SORTING LM BAND FOR '+ item)
        print('SORTING N BAND FOR '+ item)


        for file in OIFITS_TOT_N:

            AT_index = np.where(OIFITS_TOT_N[i]['FLUX']['AT_NUMBER']==item)[0][0]

            cond = np.logical_or(file['FLUX']['FLUX'][AT_index]>=FLUX_N_median[j]+2*FLUX_N_std[j],\
                                  file['FLUX']['FLUX'][AT_index]<=FLUX_N_median[j]-2*FLUX_N_std[j])


            file['FLUX']['FLAG'][AT_index] = cond | file['FLUX']['FLAG'][AT_index]
                
    
        print('END SORTING N BAND FOR '+ item)
        
        j=j+1

    

    return OIFITS_TOT_LM, OIFITS_TOT_N
