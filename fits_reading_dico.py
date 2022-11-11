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
from rhapsody_init import PROCESS_DIR, ntelescope, OIFITS_FLUX, DATA_band_name, REFLAGGING_DATA, CONCATENATE
from data_flag import FLAG_DATA, FLAG_DATA_CONCATENATE

def OIFITS_SORTING():
    
    """
        
    This routine has no direct inputs and outputs but is based on the user inputs directories and goes on the appropriate folders where the data are located.
    Here are the main steps of the routine:
        
    1) FLAG THE DATA: first the routine will look at the visibilities and closure phase in order to update the FLAG associated to visibility and closure phase
    measurements. 
    the routine will
    
    2) Then the routine will create subdirectories in the folder where the data are located and will move the fits file in the corresponding subdirectories.
       PROCESS_DIR_LM for the LM-band data, PROCESS_DIR_N for the N-band data and PROCESS_DIR_TRASH for the bad fits file according to the FLAGs.
       
    3) Then, a second sorting based on the flux will be executed to keep only the best quality data set.
     
    """    
    if CONCATENATE == True:
        FLAG_DATA_CONCATENATE()
        
    else:    
        FLAG_DATA() # If REFLAGGING_DATA == False it will only FLAG the data according to the wavelength max and min delimitation entered by the user in the init file.

    # PROCESS_DIR_TRASH   = PROCESS_DIR + '/TRASH'

    # for i in range(len(DATA_band_name)):        
        
    #     for filenames in glob.glob(PROCESS_DIR+'/'+DATA_band_name[i]+'/*.fits'):
    #         with fits.open(filenames, memmap=False) as fichier:
    #             # In the case where we are supposed to have OIFITS_FLUX table we will operate a sorting with respect to the presence or not of such table. 

    #             V2Flag_fichier= fichier["OI_VIS2"].data["FLAG"]  # ==> CHECK THIS                                       

    #             if OIFITS_FLUX == True:
                
    #                 try :
    #                     FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
    #                     fichier.close()
    #                     if np.logical_or(np.where(V2Flag_fichier == False)[0].size == 0, np.where(FLUX_MATISSE_data[:,np.where(V2Flag_fichier[0]==False)]<0)[0].size !=0):
    #                         shutil.move(filenames,PROCESS_DIR_TRASH+'/')       

    
    #                 except:
    #                     print("WARNING: %s does not include FLUX Table while it should be. This is corralted flux and not visibilities, file is moved in the TRASH folder."%filenames)
    #                     fichier.close()
    #                     shutil.move(filenames,PROCESS_DIR_TRASH+'/')       
                        
    #             else:
    #                 fichier.close()
    #                 if np.where(V2Flag_fichier == False)[0].size == 0:
    #                     shutil.move(filenames,PROCESS_DIR_TRASH+'/')       
                    
    
    #             # In the case where one of the conditions above occured, we move the file directly in the folder TRASH
                    
                                                
    return


def OIFITS_READING():
           
    """
        
    This routine will read and stock only sorted files located folder generated by the OIFITS_SORTING routine. In order to provide 
    usefull observations for the user.
    
    Outputs: O
    
    OIFITS_TOT:  Dictionnary of all the main observations available with the same keywords than in the .fits file
     
    """    

    telescope = ntelescope                

    list_files=np.zeros(len(DATA_band_name), dtype=object)
    for i in range(len(DATA_band_name)):
        PROCESS_DIR_tmp = PROCESS_DIR + '/' + DATA_band_name[i]
        list_files[i] = glob.glob(PROCESS_DIR_tmp+'/*.fits')

    OIFITS_TOT = [np.zeros(len(list_files[i]), dtype=object) for i in range(len(DATA_band_name))]

# TAKE ALL DATA AT ONCE    
    for i in range(len(list_files)):
        
            for k in range(len(list_files[i])):

                with fits.open(list_files[i][k], memmap=False) as fichier:
                    
                    # VISIBILITY
                    
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
            
                    baseline  = int(telescope*(telescope-1)/2)
                    exposure =  int(np.shape(fichier['OI_VIS2'].data['VIS2DATA'])[0]/baseline)        
                    deg = rad_to_deg(np.arctan2(u_coord_fichier,v_coord_fichier))
                    flag = flag_PA_new(deg)
                    
                    # DICTIONNAIRE
                    
                    dic                       = {'WAVEL': lambda_fichier}
                    dic['NAME']               = list_files[i][k]         
                    
                    dic['VIS2'] = {}
                    dic['VIS2']['WAVEL']      = np.array([lambda_fichier]*(int(baseline*exposure)))               
                    dic['VIS2']['BASELINE']   = np.array([baseline_fichier]*len(lambda_fichier)).T
                    dic['VIS2']['PA']         = flag
                    dic['VIS2']['U']          = np.array([u_coord_fichier]*len(lambda_fichier)).T      
                    dic['VIS2']['V']          = np.array([v_coord_fichier]*len(lambda_fichier)).T      
        
                    dic['VIS2']['VIS2']       = Visibility_2_fichier   
                    dic['VIS2']['VIS2_ERR']   = Visibility_2_err_fichier
                    dic['VIS2']['FLAG']       = V2Flag_fichier
                    
        
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

                    if OIFITS_FLUX == True:

                        FLUX_MATISSE_data = fichier['OI_FLUX'].data['FLUXDATA']
                        FLUX_MATISSE_err = fichier['OI_FLUX'].data['FLUXERR']
                        FLUX_MATISSE_flag = fichier['OI_FLUX'].data['FLAG']
                        FLUX_STA_INDEX    = fichier['OI_FLUX'].data['STA_INDEX']                
        
            
                        TEL_NAME    = fichier['OI_ARRAY'].data['TEL_NAME']                
                        STA_INDEX  = fichier['OI_ARRAY'].data['STA_INDEX']                
            
                        # GOOD_FLUX_STATION_NAME_FLUX 
                        flux_name_STA = {item: idx for idx, item in enumerate(FLUX_STA_INDEX)}
                        flux_number = ((np.array(TEL_NAME)[[flux_name_STA.get(item) for item in STA_INDEX]].astype(list)*(np.ones((len(lambda_fichier),len(FLUX_STA_INDEX)))).astype(int)).T)
                        
                        dic['FLUX'] = {}
                        dic['FLUX']['WAVEL']      = np.array([lambda_fichier]*(int(telescope)))               
                        dic['FLUX']['FLUX']       = FLUX_MATISSE_data   
                        dic['FLUX']['FLUX_ERR']   = FLUX_MATISSE_err
                        dic['FLUX']['FLAG']       = FLUX_MATISSE_flag
                        dic['FLUX']['AT_NUMBER']  = flux_number
        
                    OIFITS_TOT[i][k]   = dic
    
    
                    fichier.close()


    print('OIFITS READING OVER')
    return OIFITS_TOT


def OIFITS_READING_concatenate():
    
    list_files=np.zeros(len(DATA_band_name), dtype=object)
    for i in range(len(DATA_band_name)):
        PROCESS_DIR_tmp = PROCESS_DIR + '/' + DATA_band_name[i]
        list_files[i] = glob.glob(PROCESS_DIR_tmp+'/*.fits')

    OIFITS_TOT = [np.zeros(len(list_files[i]), dtype=object) for i in range(len(DATA_band_name))]
    
    for f in range(len(list_files)):
    
        for k in range(len(list_files[f])):
            
            with fits.open(list_files[i][k], memmap=False) as fichier:

    
                dic                       = {}
                dic['NAME']               = list_files[f]         
        
                name_HDU = np.array([fichier[t].name for t in range(len(fichier))])
                
                index_vis = np.where(name_HDU=='OI_VIS2')[0]
                
                index_wvl = np.where(name_HDU=='OI_WAVELENGTH')[0]
                
                wavel_shape = np.array([len(fichier[t].data['EFF_WAVE']) for t in index_wvl])
                
                vis2     = []
                vis2_err = []
                u_coord = []
                v_coord = []
                vis2_flag = []
                wavel     = []
                flag = []
                
                for i in range(len(index_vis)):
                    
                    if len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0])) == max(wavel_shape):
                        
                        vis2.extend(np.array(fichier[index_vis[i]].data['VIS2DATA']))
                        vis2_err.extend(np.array(fichier[index_vis[i]].data['VIS2ERR']))
                        vis2_flag.extend((np.array(fichier[index_vis[i]].data['FLAG'])))
                        u_coord.extend(np.array([fichier[index_vis[i]].data['UCOORD']]*len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0]))).T)
                        v_coord.extend(np.array([fichier[index_vis[i]].data['VCOORD']]*len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0]))).T)

                        deg = rad_to_deg(np.arctan2(np.array(fichier[index_vis[i]].data['UCOORD']),np.array(fichier[index_vis[i]].data['VCOORD'])))
                        flag.extend(flag_PA_new(deg))
                        
    
                        
                        index=index_wvl[wavel_shape==len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0]))][0] 
                        
                        wavel.extend(np.array([fichier[index].data['EFF_WAVE']]*len(fichier[index_vis[i]].data['VIS2ERR'])))
                    
                    else:
                        diff = max(wavel_shape)-len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0]))
                
                        shape = np.shape(np.array(fichier[index_vis[i]].data['VIS2DATA']))[0]
                        add = np.empty((shape,diff))*np.nan
                
                        vis2.extend(np.hstack((np.array(fichier[index_vis[i]].data['VIS2DATA']),add)))
                        vis2_err.extend(np.hstack((np.array(fichier[index_vis[i]].data['VIS2ERR']),add)))
                        vis2_flag.extend(np.hstack(((np.array(fichier[index_vis[i]].data['FLAG'])),add)))
                        u_coord.extend(np.hstack((np.array([fichier[index_vis[i]].data['UCOORD']]*len(np.array(fichier[index_vis[i]].data['VIS2DATA'])[0])).T,add)))
                        v_coord.extend(np.hstack((np.array([fichier[index_vis[i]].data['VCOORD']]*len(np.array(fichier[index_vis[i]].data['VIS2DATA'])[0])).T,add)))

                        u_temp = np.array(fichier[index_vis[i]].data['UCOORD'])
                        v_temp = np.array(fichier[index_vis[i]].data['VCOORD'])
                        deg = rad_to_deg(np.arctan2(u_temp,v_temp))
                        flag.extend(flag_PA_new(deg))
                        
                        index=index_wvl[wavel_shape==len((np.array(fichier[index_vis[i]].data['VIS2DATA'])[0]))][0] 
                        
                        wavel.extend(np.hstack((np.array([fichier[index].data['EFF_WAVE']]*len(fichier[index_vis[i]].data['VIS2ERR'])),add)))
                
                
                
                vis2     = np.array(vis2, dtype='object').astype(float)
                vis2_err = np.array(vis2_err, dtype='object').astype(float)
                vis2_flag = np.array(vis2_flag, dtype='bool') # ==> MAGIC TRICK! Like this all the NaN value that has been added to fit the wavelengths will automatically be flagged! But we will conserve the exact shape for all data, MOUHAHAHA (yes it's 2am and I'm proud of it)
                u_coord  = np.array(u_coord, dtype='object').astype(float)
                v_coord  = np.array(v_coord, dtype='object').astype(float)
                wavel    = np.array(wavel, dtype='object').astype(float)
                flag     = np.array(flag, dtype='object').astype(float)
            
                dic['VIS2'] = {}
                dic['VIS2']['WAVEL']      = wavel               
                dic['VIS2']['BASELINE']   = (u_coord**2+v_coord**2)**(1/2)
        

        
                
                dic['VIS2']['PA']         = flag
                dic['VIS2']['U']          = u_coord    
                dic['VIS2']['V']          = v_coord  
        
        
                dic['VIS2']['VIS2']       = vis2  
                dic['VIS2']['VIS2_ERR']   = vis2_err
                dic['VIS2']['FLAG']       = vis2_flag
    
                OIFITS_TOT[f][k]   = dic
    
                fichier.close()

    
    print('OIFITS READING OVER')

    return OIFITS_TOT