# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:34:08 2020

@author: jdrevon
"""
import numpy.ma as ma
from astropy.io import fits 
import sys 
import numpy as np
import matplotlib.pyplot as plt
from rhapsody_init import ntelescope

def FLUX_FLAG(OIFITS_TOT_LM, OIFITS_TOT_N, PLOT):
    
    if ntelescope == 2:

        AT=['AT1','AT2']
        
    if ntelescope == 3:

        AT=['AT1','AT2','AT3']
        
    if ntelescope == 4:
    
        AT=['AT1','AT2','AT3','AT4']
    
    for item in AT:
        
        FLUX_LM      = []
        FLUX_N       = []
        
        wavel_LM      = []
        wavel_N       = []
        
        for i in range(len(OIFITS_TOT_LM)):
            
            cond = OIFITS_TOT_LM[i]['FLUX']['AT_NUMBER']==item
            FLUX_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['FLUX'].astype('float'),mask=~cond))
            wavel_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['WAVEL'].astype('float'),mask=~cond))
                
        for i in range(len(OIFITS_TOT_N)):
            
            cond = OIFITS_TOT_N[i]['FLUX']['AT_NUMBER']==item
            FLUX_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['FLUX'].astype('float'),mask=~cond))
            wavel_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['WAVEL'].astype('float'),mask=~cond))
            
        
        FLUX_N_median  = ma.median(ma.array(FLUX_N),axis=0)
        FLUX_LM_median = ma.median(ma.array(FLUX_LM),axis=0)
    
        FLUX_N_std  = ma.median(ma.abs(ma.array(FLUX_N)-FLUX_N_median),axis=0)*1.48
        FLUX_LM_std = ma.median(ma.abs(ma.array(FLUX_LM)-FLUX_LM_median),axis=0)*1.48
        
        if PLOT==True:
        
            plt.figure()
            for i in range(len(OIFITS_TOT_LM)):
                plt.scatter(wavel_LM[i],FLUX_LM[i], c='black', s=2, alpha=0.1)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median, label = 'Mean', s=4)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median + 2*FLUX_LM_std, label = '+sigma', s=4)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median - 2*FLUX_LM_std, label = '-sigma', s=4)
            ax=plt.gca()
            ax.set_ylim([-100,5000])
            ax.set_title('LM '+ item + ' before sorting')
            plt.legend()

            plt.figure()
            for i in range(len(OIFITS_TOT_N)):
                plt.scatter(wavel_N[i],FLUX_N[i], c='black', s=2, alpha=0.1)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median, label = 'Mean', s=4)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median + 2*FLUX_N_std, label = '+sigma', s=4)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median - 2*FLUX_N_std, label = '-sigma', s=4)
            ax=plt.gca()
            ax.set_ylim([-100,1000])
            ax.set_title('N '+ item + ' before sorting')
            plt.legend()


        print('SORTING LM BAND FOR '+ item)
        for file in OIFITS_TOT_LM:
            cond = np.logical_or(file['FLUX']['FLUX'][file['FLUX']['AT_NUMBER']==item]>=FLUX_LM_median+2*FLUX_LM_std,\
                                  file['FLUX']['FLUX'][file['FLUX']['AT_NUMBER']==item]<=FLUX_LM_median-2*FLUX_LM_std)

            if np.shape(np.where(cond==True)[0])[0]!=0:
                # print(file['NAME'])
                with fits.open(file['NAME'], memmap=False, mode='update') as fichier:
    
                    TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
                    STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                
                    FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                
                    lambda_fichier       = fichier["OI_WAVELENGTH"].data["EFF_WAVE"]

                    # print(STA_INDEX)
                    # print(FLUX_STA_INDEX)
                    # print(TEL_NAME)

                    right_order = np.array([np.where(np.array(STA_INDEX)[index]== np.array(FLUX_STA_INDEX))[0][0] for index in range(len(FLUX_STA_INDEX))])
                    flux_number = ((TEL_NAME[right_order]*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)

                    # print(TEL_NAME[right_order])
                    # print(right_order)

    
                    # fichier['OI_FLUX'].data['FLAG'][flux_number==item] =\
                    # np.full(np.shape(fichier['OI_FLUX'].data['FLAG'][flux_number==item]), False)
        
    
                
                    fichier['OI_FLUX'].data['FLAG'][flux_number==item] =\
                    np.full(np.shape(fichier['OI_FLUX'].data['FLAG'][flux_number==item]), True)
    
                    file['FLUX']['FLAG'][file['FLUX']['AT_NUMBER']==item] =\
                        np.full(np.shape(file['FLUX']['FLAG'][file['FLUX']['AT_NUMBER']==item]), True)
                
                    fichier.flush()
                    fichier.close() # Au cas où mais pas nécessaire juste double sécurité de fermeture OIFITS
            

        print('END SORTING LM BAND FOR '+ item)
        print('SORTING N BAND FOR '+ item)


        for file in OIFITS_TOT_N:
            cond = np.logical_or(file['FLUX']['FLUX'][file['FLUX']['AT_NUMBER']==item]>=FLUX_N_median+2*FLUX_N_std,\
                                  file['FLUX']['FLUX'][file['FLUX']['AT_NUMBER']==item]<=FLUX_N_median-2*FLUX_N_std)

            if np.shape(np.where(cond==True)[0])[0]!=0:                                              
                with fits.open(file['NAME'], memmap=False, mode='update') as fichier:
    
                    TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
                    STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                
                    FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                
                    lambda_fichier       = fichier["OI_WAVELENGTH"].data["EFF_WAVE"]

                    right_order = np.array([np.where(np.array(STA_INDEX)[index]== np.array(FLUX_STA_INDEX))[0][0] for index in range(len(FLUX_STA_INDEX))])
                    flux_number = ((TEL_NAME[right_order]*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)
    
                    # fichier['OI_FLUX'].data['FLAG'][flux_number==item] =\
                    # np.full(np.shape(fichier['OI_FLUX'].data['FLAG'][flux_number==item]), False)
                    
                    fichier['OI_FLUX'].data['FLAG'][flux_number==item] =\
                    np.full(np.shape(fichier['OI_FLUX'].data['FLAG'][flux_number==item]), True)

                    file['FLUX']['FLAG'][file['FLUX']['AT_NUMBER']==item] =\
                        np.full(np.shape(file['FLUX']['FLAG'][file['FLUX']['AT_NUMBER']==item]), True)
            
                    fichier.flush()
                    fichier.close() # Au cas où mais pas nécessaire juste double sécurité de fermeture OIFITS
            
    
        print('END SORTING N BAND FOR '+ item)


                
        if PLOT==True:

            FLUX_LM      = []
            FLUX_N       = []
            
            wavel_LM      = []
            wavel_N       = []
            
            for i in range(len(OIFITS_TOT_LM)):
                
                cond = np.logical_and(OIFITS_TOT_LM[i]['FLUX']['FLAG']==False,OIFITS_TOT_LM[i]['FLUX']['AT_NUMBER']==item)
                FLUX_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['FLUX'].astype('float'),mask=~cond))
                wavel_LM.append(ma.masked_array(OIFITS_TOT_LM[i]['FLUX']['WAVEL'].astype('float'),mask=~cond))
                    
            for i in range(len(OIFITS_TOT_N)):
                
                cond = np.logical_and(OIFITS_TOT_N[i]['FLUX']['FLAG']==False,OIFITS_TOT_N[i]['FLUX']['AT_NUMBER']==item)
                FLUX_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['FLUX'].astype('float'),mask=~cond))
                wavel_N.append(ma.masked_array(OIFITS_TOT_N[i]['FLUX']['WAVEL'].astype('float'),mask=~cond))
            

            plt.figure()
            for i in range(len(OIFITS_TOT_LM)):
                plt.scatter(wavel_LM[i],FLUX_LM[i], c='black', s=2, alpha=0.1)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median, label = 'Mean', s=4)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median + 2*FLUX_LM_std, label = '+2sigma', s=4)
            plt.scatter(OIFITS_TOT_LM[0]['FLUX']['WAVEL'],FLUX_LM_median - 2*FLUX_LM_std, label = '-2sigma', s=4)
            ax=plt.gca()
            ax.set_ylim([-100,5000])
            ax.set_title('LM '+ item + ' after sorting')
            plt.legend()

            plt.figure()
            for i in range(len(OIFITS_TOT_N)):
                plt.scatter(wavel_N[i],FLUX_N[i], c='black', s=2, alpha=0.1)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median, label = 'Mean', s=4)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median + 2*FLUX_N_std, label = '+2sigma', s=4)
            plt.scatter(OIFITS_TOT_N[0]['FLUX']['WAVEL'],FLUX_N_median - 2*FLUX_N_std, label = '-2sigma', s=4)
            ax=plt.gca()
            ax.set_ylim([-100,1000])
            ax.set_title('N '+ item + ' after sorting')
            plt.legend()



    return
