#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:34:08 2020

@author: jdrevon
"""
import numpy as np
from rhapsody_init import ntelescope


def FLUX_FLAG(OIFITS_TOT):
    
    if ntelescope == 2:

        AT=['AT1','AT2']
        
    if ntelescope == 3:

        AT=['AT1','AT2','AT3']
        
    if ntelescope == 4:
    
        AT=['AT1','AT2','AT3','AT4']
    
        
    # What we want is to make statistics per AT on the flux to flag all the flux data
    # So we will stock all the data per AT make the median of all the flux in order to have N median spectras for the N AT's 
    
    
    # I initiate the array which will contain the information for each band, for each ATs, the flux and the wavelengths values associated.    


    for k in range(len(OIFITS_TOT)): # Loop over all the bandwidth 

        FLUX_stock = [np.zeros((len(OIFITS_TOT[k]),ntelescope,len(OIFITS_TOT[k][0]['WAVEL'])), dtype=object) for i in range(len(OIFITS_TOT))]
        WAVEL_stock = [np.zeros((len(OIFITS_TOT[k]),ntelescope,len(OIFITS_TOT[k][0]['WAVEL'])), dtype=object) for i in range(len(OIFITS_TOT))]
    
        for j in range(len(AT)): # Loop over the ATs of the k-th bandwidth
    
            for i in range(len(OIFITS_TOT[k])): # Loop over the files of the k-th bandwidth 
                
                # I select the flux data from the j-th AT of the k-th BAND in the i-th file
                
                AT_index = np.where(OIFITS_TOT[k][i]['FLUX']['AT_NUMBER']==AT[j])[0][0]
                
                # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
                # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
                            
                FLUX_stock[k][i][j]  = OIFITS_TOT[k][i]['FLUX']['FLUX'][AT_index]
                
                WAVEL_stock[k][i][j] = OIFITS_TOT[k][i]['FLUX']['WAVEL'][AT_index]

                

            # Now it's time to start to make statistics on the flux for each band data that we have stocked  for a given AT
            # We take for each band the median value at a given AT and compute the standard deviation
            # Then we can compare each flux value for each file of each bands and reflag
            
            # I want to make the median over the flux from the k-th band with the j-th telescope
            
            FLUX_median  = np.median(np.array(FLUX_stock[k]),axis=1)
            FLUX_std     = np.median(np.abs(np.array(FLUX_stock[k])-FLUX_median[k]),axis=1)*1.48

            # Now we have to loop over each file one more time to change the FLAG accordingly to the following conditions written in the variable cond            

            for file in OIFITS_TOT[k]:
    
                AT_index = np.where(file['FLUX']['AT_NUMBER']==AT[j])[0][0]
    
                cond = np.logical_or(file['FLUX']['FLUX'][AT_index]>=FLUX_median[j]+2*FLUX_std[j],\
                                      file['FLUX']['FLUX'][AT_index]<=FLUX_median[j]-2*FLUX_std[j])
    
    
                file['FLUX']['FLAG'][AT_index] = cond | file['FLUX']['FLAG'][AT_index]
                    
    
            print('END SORTING BAND FOR '+ AT[j])


    

    return OIFITS_TOT



def FLUX_FLAG_CONCATENATE(OIFITS_TOT):
    
    if ntelescope == 2:

        AT=['AT1','AT2']
        
    if ntelescope == 3:

        AT=['AT1','AT2','AT3']
        
    if ntelescope == 4:
    
        AT=['AT1','AT2','AT3','AT4']
    
        
    # What we want is to make statistics per AT on the flux to flag all the flux data
    # So we will stock all the data per AT make the median of all the flux in order to have N median spectras for the N AT's 
    
    
    # I initiate the array which will contain the information for each band, for each ATs, the flux and the wavelengths values associated.    


    FLUX_stock  = np.zeros((ntelescope, int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope), len(OIFITS_TOT[0]['WAVEL'])), dtype=object)
    WAVEL_stock = np.zeros((ntelescope, int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope), len(OIFITS_TOT[0]['WAVEL'])), dtype=object)             
    
    for j in range(len(AT)): # Loop over the ATs of the k-th bandwidth

        
        # I select the flux data from the j-th AT of the k-th BAND in the i-th file
        
        AT_index = np.where(OIFITS_TOT[0]['FLUX']['AT_NUMBER']==AT[j])
        
        # Once I have the correspond index of the flux from AT over the N arrays of flux I will stock the values in FLUX_LM[i][j][:]
        # So in the first loop AT1 values will be stocked in the first slot of the FLUX_LM array.
                    
        FLUX_stock[j]  = np.reshape(OIFITS_TOT[0]['FLUX']['FLUX'][AT_index], (int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope),len(OIFITS_TOT[0]['WAVEL'])))        
        WAVEL_stock[j] = np.reshape(OIFITS_TOT[0]['FLUX']['WAVEL'][AT_index], (int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope),len(OIFITS_TOT[0]['WAVEL'])))
            

        # Now it's time to start to make statistics on the flux for each band data that we have stocked  for a given AT
        # We take for each band the median value at a given AT and compute the standard deviation
        # Then we can compare each flux value for each file of each bands and reflag
        
        
    FLUX_median  = np.median(np.array(FLUX_stock),axis=1)
    FLUX_std     = [np.median(np.abs(np.array(FLUX_stock[j])-FLUX_median[j]),axis=0)*1.48 for i in range(ntelescope)]
    
    for j in range(len(AT)): # Loop over the ATs of the k-th bandwidth

        file = OIFITS_TOT[0]
        AT_index = np.where(file['FLUX']['AT_NUMBER']==AT[j])
    
        FLUX_M_mod = [FLUX_median[j]]*int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope)
        FLUX_M_mod = np.reshape(FLUX_M_mod,(np.shape(FLUX_M_mod)[0]*np.shape(FLUX_M_mod)[1]))
    
        FLUX_std_mod = [FLUX_std[j]]*int(len(OIFITS_TOT[0]['FLUX']["FLUX"])/ntelescope)
        FLUX_std_mod = np.reshape(FLUX_std_mod,(np.shape(FLUX_std_mod)[0]*np.shape(FLUX_std_mod)[1]))
    
    
        cond = np.logical_or(file['FLUX']['FLUX'][AT_index]>=FLUX_M_mod[j]+2*FLUX_std_mod[j],\
                              file['FLUX']['FLUX'][AT_index]<=FLUX_M_mod[j]-2*FLUX_std_mod[j])
    
    
        file['FLUX']['FLAG'][AT_index] = cond | file['FLUX']['FLAG'][AT_index]
                
    
        print('END SORTING BAND FOR '+ AT[j])


    

    return OIFITS_TOT

