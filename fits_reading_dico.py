#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 14:21:51 2020

@author: jdrevon
"""
import numpy as np
import os 
from astropy.io import fits 
import glob
from A7_rad_to_deg import rad_to_deg
from A8_PA_flag_new import flag_PA_new
import shutil
from rhapsody_init import MATISSE_DIR, ntelescope
from data_flag import FLAG_DATA
from data_flag_flux import FLUX_FLAG

def OIFITS_SORTING():
    
    """
        
    This routine has no direct inputs and outputs but is based on the user inputs directories and goes on the appropriate folders where the data are located.
    Here are the main steps of the routine:
        
    1) FLAG THE DATA: first the routine will look at the visibilities and closure phase in order to update the FLAG associated to visibility and closure phase
    measurements. 
    the routine will
    
    2) Then the routine will create subdirectories in the folder where the data are located and will move the fits file in the corresponding subdirectories.
       MATISSE_DIR_LM for the LM-band data, MATISSE_DIR_N for the N-band data and MATISSE_DIR_TRASH for the bad fits file according to the FLAGs.
       
    3) Then, a second sorting based on the flux will be executed to keep only the best quality data set.
     
    """    
    
    print('START FLAGGING')

    FLAG_DATA()

    print('CREATING NEW DIRECTORIES')

    MATISSE_DIR_LM      = MATISSE_DIR + '/LM' 
    MATISSE_DIR_N       = MATISSE_DIR + '/N'
    MATISSE_DIR_TRASH   = MATISSE_DIR + '/TRASH'


    # move file in correct folder
    if not os.path.exists(MATISSE_DIR_LM):
        os.makedirs(MATISSE_DIR_LM)
    if not os.path.exists(MATISSE_DIR_N):
        os.makedirs(MATISSE_DIR_N)

    if not os.path.exists(MATISSE_DIR_TRASH):
        os.makedirs(MATISSE_DIR_TRASH)

    
    for filenames in glob.glob(MATISSE_DIR+'/*.fits'):

        with fits.open(filenames, memmap=False) as fichier:
            FLAG = False
            
            try: 
                FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
                V2Flag_fichier= fichier["OI_VIS2"].data["FLAG"]
                if np.logical_or(np.where(V2Flag_fichier == False)[0].size == 0, np.where(FLUX_MATISSE_data<0)[0].size !=0):
                    FLAG = True
                    
            except : 
                FLAG = True

            fichier.close()
            if FLAG == True:
                
                shutil.move(filenames,MATISSE_DIR_TRASH+'/')                
                            
            elif np.logical_and(FLAG == False,'IR-N' in filenames):
                shutil.move(filenames,MATISSE_DIR_N+'/')
            
            elif np.logical_and(FLAG == False,'IR-LM' in filenames):                   
                shutil.move(filenames,MATISSE_DIR_LM+'/')


    list_files_LM = (glob.glob(MATISSE_DIR_LM+'/*.fits'))
    list_files_N = (glob.glob(MATISSE_DIR_N+'/*.fits'))
    count_N = 0
    count_LM = 0
    OIFITS_TOT_LM = np.zeros(len(list_files_LM), dtype = 'object')
    OIFITS_TOT_N = np.zeros(len(list_files_N), dtype = 'object')

    list_files = (glob.glob(MATISSE_DIR_LM+'/*.fits')+ glob.glob(MATISSE_DIR_N+'/*.fits'))

    print('START READING DATA FOR FLUX SORTING')
    for filenames in list_files:

        with fits.open(filenames, memmap=False) as fichier:

            
            Visibility_2_fichier = fichier["OI_VIS2"].data["VIS2DATA"]
            Visibility_2_err_fichier = fichier["OI_VIS2"].data["VIS2ERR"]
            V2Flag_fichier       = fichier["OI_VIS2"].data["FLAG"]
            u_coord_fichier      = fichier["OI_VIS2"].data["UCOORD"]
            v_coord_fichier      = fichier["OI_VIS2"].data["VCOORD"]
            baseline_fichier     = np.sqrt(np.array(u_coord_fichier)**2+np.array(v_coord_fichier)**2).tolist()
            lambda_fichier       = fichier["OI_WAVELENGTH"].data["EFF_WAVE"]


            # CLOSURE PHASE            
            
            U1 = fichier["OI_T3"].data["U1COORD"]
            U2 = fichier["OI_T3"].data["U2COORD"]
            V1 = fichier["OI_T3"].data["V1COORD"]
            V2 = fichier["OI_T3"].data["V2COORD"]
            U3 = U1+U2
            V3 = V1+V2
            T3_PHI = fichier["OI_T3"].data["T3PHI"]
            T3_PHI_ERR = fichier["OI_T3"].data["T3PHIERR"]
            FLAG_T3 = fichier["OI_T3"].data["FLAG"]
    
            telescope = ntelescope
            baseline  = int(telescope*(telescope-1)/2)
            exposure =  int(np.shape(fichier['OI_VIS2'].data['VIS2DATA'])[0]/baseline)        

            FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
            FLUX_MATISSE_err = fichier['OI_FLUX'].data['FLUXERR']
            FLUX_MATISSE_flag = fichier['OI_FLUX'].data['FLAG']
            FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                

            TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
            STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                

            # GOOD_FLUX_STATION_NAME_FLUX 

            flux_name_STA = {item: idx for idx, item in enumerate(FLUX_STA_INDEX)}
            flux_number = ((np.array(TEL_NAME)[[flux_name_STA.get(item) for item in STA_INDEX]].astype(list)*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)   
            deg = rad_to_deg(np.arctan2(u_coord_fichier,v_coord_fichier))
            flag = flag_PA_new(deg)

            # DICTIONNAIRE
            
            dic                       = {'WAVEL': lambda_fichier}
            dic['NAME']               = filenames         
            
            dic['VIS2'] = {}
            dic['VIS2']['WAVEL']      = np.array([lambda_fichier]*(int(baseline*exposure)))               
            dic['VIS2']['BASELINE']   = np.array([baseline_fichier]*len(lambda_fichier)).T
            dic['VIS2']['U']          = u_coord_fichier 
            dic['VIS2']['V']          = v_coord_fichier 
            dic['VIS2']['VIS2']       = Visibility_2_fichier   
            dic['VIS2']['VIS2_ERR']   = Visibility_2_err_fichier
            dic['VIS2']['PA']         = flag
            dic['VIS2']['FLAG']       = V2Flag_fichier
                
            dic['FLUX'] = {}
            dic['FLUX']['WAVEL']      = np.array([lambda_fichier]*(int(telescope)))               
            dic['FLUX']['FLUX']       = FLUX_MATISSE_data   
            dic['FLUX']['FLUX_ERR']   = FLUX_MATISSE_err
            dic['FLUX']['FLAG']       = FLUX_MATISSE_flag
            dic['FLUX']['AT_NUMBER']  = flux_number

            dic['T3'] = {}
            dic['T3']['WAVEL']        = np.array([lambda_fichier]*(int(telescope*exposure)))               
            dic['T3']['T3']           = T3_PHI   
            dic['T3']['T3_ERR']       = T3_PHI_ERR
            dic['T3']['B1']           = np.array([np.sqrt(U1**2+V1**2)]*int(len(lambda_fichier))).T
            dic['T3']['B2']           = np.array([np.sqrt(U2**2+V2**2)]*int(len(lambda_fichier))).T
            dic['T3']['B3']           = np.array([np.sqrt(U3**2+V3**2)]*int(len(lambda_fichier))).T
            dic['T3']['FLAG']         = FLAG_T3


            if '-N_' in filenames:
                
                OIFITS_TOT_N[count_N]   = dic
                count_N = count_N+1
                
            else:
                
                OIFITS_TOT_LM[count_LM]   = dic
                count_LM = count_LM+1
        
    print('STARTING FLUX SORTING')
    FLUX_FLAG(OIFITS_TOT_LM, OIFITS_TOT_N, False)
    
    return
    

def OIFITS_READING_complete():
           
    """
        
    This routine will read and stock only sorted files located in LM and N band MATISSE folder generated by the OIFITS_SORTING routine. In order to provide 
    usefull observations for the user.
    
    Outputs: O
    
    V2_MATISSE : 1) Wavelengths[m], 2) Baseline[m], 3) Visibilities, 4) Error on Visibilities, 5) Position Angle[°]
    UV : 1) Wavelengths[m], 2) U[m], 3) V[m] ==> UV coordinates for the visibility data points
    UV_TP : 1) U1[m], 2) U2[m], 3) U3[m], 4) V1[m], 5) V2[m], 6) V3[m], 7) Wavelengths[m] ==> UV coordinates for the triple products
    TP_MATISSE : 1) Wavelenghts[m], 2) Triple Product, 3) Triple Product associated error bars
    FLUX_MATISSE_LM : 1) Wavelengths[m], 2) Flux MATISSE, 3) Error on the Flux [Jy], 4) Telescope Number (AT1, AT2, AT3, AT4), 5) FLAG on the data flux
    FLUX_MATISSE_N :1) Wavelengths[m], 2) Flux MATISSE, 3) Error on the Flux [Jy], 4) Telescope Number (AT1, AT2, AT3, AT4), 5) FLAG on the data flux
    BAND_WIDTH_MATISSE : 1) Wavelenghts [m], Bandwidth [m]
    OIFITS_TOT_LM : Dictionnary of all the main observations available with the same keywords than in the .fits file
    OIFITS_TOT_N :  Dictionnary of all the main observations available with the same keywords than in the .fits file
     
    """    

    telescope = ntelescope

    MATISSE_DIR_LM      = MATISSE_DIR + '/LM' 
    MATISSE_DIR_N       = MATISSE_DIR + '/N'
                
    try:
        list_files_LM = (glob.glob(MATISSE_DIR_LM+'/*.fits'))
    except: 
        list_files_LM = []
    try :
        list_files_N = (glob.glob(MATISSE_DIR_N+'/*.fits'))
    except: 
        list_files_N = []

    try:
        OIFITS_TOT_LM = np.zeros(len(list_files_LM), dtype = 'object')
    except : 
        OIFITS_TOT_LM = []
    try :    

        OIFITS_TOT_N = np.zeros(len(list_files_N), dtype = 'object')
    except : 
        OIFITS_TOT_N = []

    try :    
        with fits.open(list_files_LM[0], memmap=False) as fichier:
            wavel_number_LM = fichier['OI_WAVELENGTH'].header['NAXIS2']
        
        FLUX_MATISSE_LM         = np.zeros((len(list_files_LM),5,telescope ,wavel_number_LM),dtype = 'object')
    except : 
        FLUX_MATISSE_LM         = []

    try :    
        with fits.open(list_files_N[0], memmap=False) as fichier:
            wavel_number_N = fichier['OI_WAVELENGTH'].header['NAXIS2']
        FLUX_MATISSE_N          = np.zeros((len(list_files_N),5,telescope ,wavel_number_N),dtype = 'object')
    except : 
        FLUX_MATISSE_N          = []
    

    V2_MATISSE              = np.zeros((1,5)) # [lambda,base,V[lambda,base],V_erreur[lambda,base],PA]
    UV                      = np.zeros((1,3))
    UV_TP                   = np.zeros((1,7))
    TP_MATISSE              = np.zeros((1,3))
    BAND_WIDTH_MATISSE      = np.zeros((1,2))



    list_files = (list_files_LM + list_files_N)

# TAKE ALL DATA AT ONCE    

    count_N = 0
    count_LM = 0
    
    
    for filenames in list_files:

        with fits.open(filenames, memmap=False) as fichier:
            
            # VISIBILITY
            
            Visibility_2_fichier = fichier["OI_VIS2"].data["VIS2DATA"]
            Visibility_2_err_fichier = fichier["OI_VIS2"].data["VIS2ERR"]
            V2Flag_fichier       = fichier["OI_VIS2"].data["FLAG"]
            u_coord_fichier      = fichier["OI_VIS2"].data["UCOORD"]
            v_coord_fichier      = fichier["OI_VIS2"].data["VCOORD"]
            baseline_fichier     = np.sqrt(np.array(u_coord_fichier)**2+np.array(v_coord_fichier)**2).tolist()
            lambda_fichier       = fichier["OI_WAVELENGTH"].data["EFF_WAVE"]
            bandwith_fichier     = fichier["OI_WAVELENGTH"].data["EFF_BAND"]


            # CLOSURE PHASE            
            
            U1 = fichier["OI_T3"].data["U1COORD"]
            U2 = fichier["OI_T3"].data["U2COORD"]
            V1 = fichier["OI_T3"].data["V1COORD"]
            V2 = fichier["OI_T3"].data["V2COORD"]
            U3 = U1+U2
            V3 = V1+V2
            T3_PHI = fichier["OI_T3"].data["T3PHI"]
            T3_PHI_ERR = fichier["OI_T3"].data["T3PHIERR"]
            FLAG_T3 = fichier["OI_T3"].data["FLAG"]
    
            baseline  = int(telescope*(telescope-1)/2)
            exposure =  int(np.shape(fichier['OI_VIS2'].data['VIS2DATA'])[0]/baseline)        
    
            BAND_WIDTH_MATISSE = np.append(BAND_WIDTH_MATISSE,np.array([lambda_fichier,bandwith_fichier]).T,axis=0)
            
            FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
            FLUX_MATISSE_err = fichier['OI_FLUX'].data['FLUXERR']
            FLUX_MATISSE_flag = fichier['OI_FLUX'].data['FLAG']
            FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                

            TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
            STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                

            # GOOD_FLUX_STATION_NAME_FLUX 
            flux_name_STA = {item: idx for idx, item in enumerate(FLUX_STA_INDEX)}
            flux_number = ((np.array(TEL_NAME)[[flux_name_STA.get(item) for item in STA_INDEX]].astype(list)*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)
            flux_stock = [np.array([lambda_fichier]*(int(telescope))),FLUX_MATISSE_data,FLUX_MATISSE_err,flux_number, FLUX_MATISSE_flag]
            deg = rad_to_deg(np.arctan2(u_coord_fichier,v_coord_fichier))
            flag = flag_PA_new(deg)
    

                
            # DICTIONNAIRE
            
            dic                       = {'WAVEL': lambda_fichier}
            dic['NAME']               = filenames         
            
            dic['VIS2'] = {}
            dic['VIS2']['WAVEL']      = np.array([lambda_fichier]*(int(baseline*exposure)))               
            dic['VIS2']['BASELINE']   = np.array([baseline_fichier]*len(lambda_fichier)).T
            dic['VIS2']['PA']         = flag
            dic['VIS2']['U']          = np.array([u_coord_fichier]*len(lambda_fichier)).T      
            dic['VIS2']['V']          = np.array([v_coord_fichier]*len(lambda_fichier)).T      


            dic['VIS2']['VIS2']       = Visibility_2_fichier   
            dic['VIS2']['VIS2_ERR']   = Visibility_2_err_fichier
            dic['VIS2']['FLAG']       = V2Flag_fichier
                
            dic['FLUX'] = {}
            dic['FLUX']['WAVEL']      = np.array([lambda_fichier]*(int(telescope)))               
            dic['FLUX']['FLUX']       = FLUX_MATISSE_data   
            dic['FLUX']['FLUX_ERR']   = FLUX_MATISSE_err
            dic['FLUX']['FLAG']       = FLUX_MATISSE_flag
            dic['FLUX']['AT_NUMBER']  = flux_number

            dic['T3'] = {}
            dic['T3']['WAVEL']        = np.array([lambda_fichier]*(int(telescope*exposure)))               
            dic['T3']['T3']           = T3_PHI   
            dic['T3']['T3_ERR']       = T3_PHI_ERR
            
            # print('HELLO THERE')

            dic['T3']['U1']           = np.array([U1]*int(len(lambda_fichier))).T
            dic['T3']['U2']           = np.array([U2]*int(len(lambda_fichier))).T
            dic['T3']['U3']           = np.array([U3]*int(len(lambda_fichier))).T

            # print('GENERAL KENOBI')


            dic['T3']['V1']           = np.array([V1]*int(len(lambda_fichier))).T
            dic['T3']['V2']           = np.array([V2]*int(len(lambda_fichier))).T
            dic['T3']['V3']           = np.array([V3]*int(len(lambda_fichier))).T

            dic['T3']['B1']           = np.array([np.sqrt(U1**2+V1**2)]*int(len(lambda_fichier))).T
            dic['T3']['B2']           = np.array([np.sqrt(U2**2+V2**2)]*int(len(lambda_fichier))).T
            dic['T3']['B3']           = np.array([np.sqrt(U3**2+V3**2)]*int(len(lambda_fichier))).T
            dic['T3']['FLAG']         = FLAG_T3

            if '-N_' in filenames:
                
                FLUX_MATISSE_N[count_N] = flux_stock
                OIFITS_TOT_N[count_N]   = dic
                count_N = count_N+1
                
            else:
                
                FLUX_MATISSE_LM[count_LM] = flux_stock 
                OIFITS_TOT_LM[count_LM]   = dic
                count_LM = count_LM+1

    
            phi_stock = [np.array([lambda_fichier]*(int(telescope*exposure)))[np.invert(FLAG_T3)],T3_PHI[np.invert(FLAG_T3)],T3_PHI_ERR[np.invert(FLAG_T3)]]
            phi_transpo = np.array(phi_stock).T
            phi_sorted = phi_transpo[(-phi_transpo[:,0]).argsort()]
            TP_MATISSE = np.append(TP_MATISSE,phi_sorted,axis=0)
            
            
            uv_tp_stock = [np.array([U1]*len(lambda_fichier)).T[np.invert(FLAG_T3)],np.array([U2]*len(lambda_fichier)).T[np.invert(FLAG_T3)],np.array([U3]*len(lambda_fichier)).T[np.invert(FLAG_T3)], \
                           np.array([V1]*len(lambda_fichier)).T[np.invert(FLAG_T3)],np.array([V2]*len(lambda_fichier)).T[np.invert(FLAG_T3)],np.array([V3]*len(lambda_fichier)).T[np.invert(FLAG_T3)],\
                               np.array([lambda_fichier]*int(telescope*exposure))[np.invert(FLAG_T3)]]
            
            uv_tp_transpo = np.array(uv_tp_stock).T   
            UV_TP = np.append(UV_TP,uv_tp_transpo,axis=0)
                    

            # v2_stock = [np.array([lambda_fichier]*(int(baseline*exposure)))[np.invert(V2Flag_fichier)], np.array([baseline_fichier]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)] ,\
            #             Visibility_2_fichier[np.invert(V2Flag_fichier)],Visibility_2_err_fichier[np.invert(V2Flag_fichier)], \
            #             np.array([deg]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)],np.array([flag]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)]]
            v2_stock = [np.array([lambda_fichier]*(int(baseline*exposure)))[np.invert(V2Flag_fichier)], np.array([baseline_fichier]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)] ,\
                        Visibility_2_fichier[np.invert(V2Flag_fichier)],Visibility_2_err_fichier[np.invert(V2Flag_fichier)], \
                        np.array([deg]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)]]
            v2_transpo = np.array(v2_stock).T
            V2_MATISSE = np.append(V2_MATISSE,v2_transpo,axis=0)
                
            uv_stock = [np.array([lambda_fichier]*int(baseline*exposure))[np.invert(V2Flag_fichier)], np.array([u_coord_fichier]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)],\
                        np.array([v_coord_fichier]*len(lambda_fichier)).T[np.invert(V2Flag_fichier)]]
    
            uv_transpo = np.array(uv_stock).T
            UV = np.append(UV,uv_transpo,axis=0)
    
    
            fichier.close()
            

    V2_MATISSE          = np.delete(V2_MATISSE,0,0) 
    UV                  = np.delete(UV,0,0) 
    UV_TP               = np.delete(UV_TP,0,0) 
    TP_MATISSE          = np.delete(TP_MATISSE,0,0)
    BAND_WIDTH_MATISSE  = np.delete(BAND_WIDTH_MATISSE,0,0)

    print('P0_OK')
    return V2_MATISSE, UV, UV_TP, TP_MATISSE, FLUX_MATISSE_LM,FLUX_MATISSE_N, BAND_WIDTH_MATISSE, OIFITS_TOT_LM, OIFITS_TOT_N


def OIFITS_READING():
           
    """
        
    This routine will read and stock only sorted files located in LM and N band MATISSE folder generated by the OIFITS_SORTING routine. In order to provide 
    usefull observations for the user.
    
    Outputs: O
    
    OIFITS_TOT_LM : Dictionnary of all the main observations available with the same keywords than in the .fits file
    OIFITS_TOT_N :  Dictionnary of all the main observations available with the same keywords than in the .fits file
     
    """    

    telescope = ntelescope

    MATISSE_DIR_LM      = MATISSE_DIR + '/LM' 
    MATISSE_DIR_N       = MATISSE_DIR + '/N'
                
    try:
        list_files_LM = (glob.glob(MATISSE_DIR_LM+'/*.fits'))
    except: 
        list_files_LM = []
    try :
        list_files_N = (glob.glob(MATISSE_DIR_N+'/*.fits'))
    except: 
        list_files_N = []

    try:
        OIFITS_TOT_LM = np.zeros(len(list_files_LM), dtype = 'object')
    except : 
        OIFITS_TOT_LM = []
    try :    

        OIFITS_TOT_N = np.zeros(len(list_files_N), dtype = 'object')
    except : 
        OIFITS_TOT_N = []

    try :    
        with fits.open(list_files_LM[0], memmap=False) as fichier:
            wavel_number_LM = fichier['OI_WAVELENGTH'].header['NAXIS2']
        
        FLUX_MATISSE_LM         = np.zeros((len(list_files_LM),5,telescope ,wavel_number_LM),dtype = 'object')
    except : 
        FLUX_MATISSE_LM         = []

    try :    
        with fits.open(list_files_N[0], memmap=False) as fichier:
            wavel_number_N = fichier['OI_WAVELENGTH'].header['NAXIS2']
        FLUX_MATISSE_N          = np.zeros((len(list_files_N),5,telescope ,wavel_number_N),dtype = 'object')
    except : 
        FLUX_MATISSE_N          = []


    list_files = (list_files_LM + list_files_N)

# TAKE ALL DATA AT ONCE    

    count_N = 0
    count_LM = 0
    
    
    for filenames in list_files:

        with fits.open(filenames, memmap=False) as fichier:
            
            # VISIBILITY
            
            Visibility_2_fichier = fichier["OI_VIS2"].data["VIS2DATA"]
            Visibility_2_err_fichier = fichier["OI_VIS2"].data["VIS2ERR"]
            V2Flag_fichier       = fichier["OI_VIS2"].data["FLAG"]
            u_coord_fichier      = fichier["OI_VIS2"].data["UCOORD"]
            v_coord_fichier      = fichier["OI_VIS2"].data["VCOORD"]
            baseline_fichier     = np.sqrt(np.array(u_coord_fichier)**2+np.array(v_coord_fichier)**2).tolist()
            lambda_fichier       = fichier["OI_WAVELENGTH"].data["EFF_WAVE"]
            bandwith_fichier     = fichier["OI_WAVELENGTH"].data["EFF_BAND"]


            # CLOSURE PHASE            
            
            U1 = fichier["OI_T3"].data["U1COORD"]
            U2 = fichier["OI_T3"].data["U2COORD"]
            V1 = fichier["OI_T3"].data["V1COORD"]
            V2 = fichier["OI_T3"].data["V2COORD"]
            U3 = U1+U2
            V3 = V1+V2
            T3_PHI = fichier["OI_T3"].data["T3PHI"]
            T3_PHI_ERR = fichier["OI_T3"].data["T3PHIERR"]
            FLAG_T3 = fichier["OI_T3"].data["FLAG"]
    
            baseline  = int(telescope*(telescope-1)/2)
            exposure =  int(np.shape(fichier['OI_VIS2'].data['VIS2DATA'])[0]/baseline)        
                
            FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
            FLUX_MATISSE_err = fichier['OI_FLUX'].data['FLUXERR']
            FLUX_MATISSE_flag = fichier['OI_FLUX'].data['FLAG']
            FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                

            TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
            STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                

            # GOOD_FLUX_STATION_NAME_FLUX 
            flux_name_STA = {item: idx for idx, item in enumerate(FLUX_STA_INDEX)}
            flux_number = ((np.array(TEL_NAME)[[flux_name_STA.get(item) for item in STA_INDEX]].astype(list)*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)
            flux_stock = [np.array([lambda_fichier]*(int(telescope))),FLUX_MATISSE_data,FLUX_MATISSE_err,flux_number, FLUX_MATISSE_flag]
            deg = rad_to_deg(np.arctan2(u_coord_fichier,v_coord_fichier))
            flag = flag_PA_new(deg)
                
            # DICTIONNAIRE
            
            dic                       = {'WAVEL': lambda_fichier}
            dic['NAME']               = filenames         
            
            dic['VIS2'] = {}
            dic['VIS2']['WAVEL']      = np.array([lambda_fichier]*(int(baseline*exposure)))               
            dic['VIS2']['BASELINE']   = np.array([baseline_fichier]*len(lambda_fichier)).T
            dic['VIS2']['PA']         = flag
            dic['VIS2']['U']          = np.array([u_coord_fichier]*len(lambda_fichier)).T      
            dic['VIS2']['V']          = np.array([v_coord_fichier]*len(lambda_fichier)).T      


            dic['VIS2']['VIS2']       = Visibility_2_fichier   
            dic['VIS2']['VIS2_ERR']   = Visibility_2_err_fichier
            dic['VIS2']['FLAG']       = V2Flag_fichier
                
            dic['FLUX'] = {}
            dic['FLUX']['WAVEL']      = np.array([lambda_fichier]*(int(telescope)))               
            dic['FLUX']['FLUX']       = FLUX_MATISSE_data   
            dic['FLUX']['FLUX_ERR']   = FLUX_MATISSE_err
            dic['FLUX']['FLAG']       = FLUX_MATISSE_flag
            dic['FLUX']['AT_NUMBER']  = flux_number

            dic['T3'] = {}
            dic['T3']['WAVEL']        = np.array([lambda_fichier]*(int(telescope*exposure)))               
            dic['T3']['T3']           = T3_PHI   
            dic['T3']['T3_ERR']       = T3_PHI_ERR
            
            # print('HELLO THERE')

            dic['T3']['U1']           = np.array([U1]*int(len(lambda_fichier))).T
            dic['T3']['U2']           = np.array([U2]*int(len(lambda_fichier))).T
            dic['T3']['U3']           = np.array([U3]*int(len(lambda_fichier))).T

            # print('GENERAL KENOBI')


            dic['T3']['V1']           = np.array([V1]*int(len(lambda_fichier))).T
            dic['T3']['V2']           = np.array([V2]*int(len(lambda_fichier))).T
            dic['T3']['V3']           = np.array([V3]*int(len(lambda_fichier))).T

            dic['T3']['B1']           = np.array([np.sqrt(U1**2+V1**2)]*int(len(lambda_fichier))).T
            dic['T3']['B2']           = np.array([np.sqrt(U2**2+V2**2)]*int(len(lambda_fichier))).T
            dic['T3']['B3']           = np.array([np.sqrt(U3**2+V3**2)]*int(len(lambda_fichier))).T
            dic['T3']['FLAG']         = FLAG_T3

            if '-N_' in filenames:
                
                FLUX_MATISSE_N[count_N] = flux_stock
                OIFITS_TOT_N[count_N]   = dic
                count_N = count_N+1
                
            else:
                
                FLUX_MATISSE_LM[count_LM] = flux_stock 
                OIFITS_TOT_LM[count_LM]   = dic
                count_LM = count_LM+1
    
    
            fichier.close()


    print('OIFITS READING OVER')
    return OIFITS_TOT_LM, OIFITS_TOT_N

