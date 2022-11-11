#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:52:31 2020

@author: jdrevon
"""

# Import stuff
from astropy.io import fits
import numpy as np
import glob
from rhapsody_init import PROCESS_DIR, DATA_band_name, DATA_band_min, DATA_band_max, OIFITS_FLUX, REFLAGGING_DATA
import sys

def FLAG_DATA():

    for k in range(len(DATA_band_name)):
    
        for filenames in glob.glob(PROCESS_DIR+'/'+DATA_band_name[k]+'/*.fits'):
                    
            wlmin = DATA_band_min[k]*1E-6
            wlmax = DATA_band_max[k]*1E-6
                    
        
            with fits.open(filenames, memmap=False, mode='update') as fichier:
    
                # print(filenames)
                # name = np.zeros(len(data), dtype=object)
                # for i in range(len(data)):
                #     name[i] = data[i].name
    
                # data=fits.open(filenames, mode='update', memmap=False)
    
                # Read data
                VIS2     = fichier['OI_VIS2'].data['VIS2DATA']
                VIS2ERR  = fichier['OI_VIS2'].data['VIS2ERR']
                
                isnan = VIS2ERR != VIS2ERR
                
                VIS2ERR[isnan] = 1E99
                
                WLEN = np.array([fichier['OI_WAVELENGTH'].data['EFF_WAVE']]*np.shape(VIS2)[0])
                # VIS2FLAG = data['OI_VIS2'].data['FLAG']
            
                vis    = np.sqrt(np.abs(VIS2)) * np.sign(VIS2)
                viserr = 0.5* VIS2ERR / (vis+(vis==0));
                
                if REFLAGGING_DATA == True:
                    
                    print('START FLAGGING')

                    
                    flag = (WLEN < wlmax)      &\
                           (WLEN > wlmin)      &\
                            (vis > 0. - viserr) &\
                            (vis < 1. + viserr) &\
                            (vis > -0.1)         &\
                            (vis < 1.1)          &\
                            (viserr > 0)         &\
                            (viserr < 0.1)    #   &\

                else:
                    
                    flag = (WLEN < wlmax)      &\
                           (WLEN > wlmin)

                        
                flag = ~flag
                try:
                    fichier['OI_VIS2'].data['FLAG'] = flag    
                except:
                    fichier['OI_VIS2'].data['FLAG'] = flag.reshape((len(flag),))
                    
                
                fichier.flush()
                fichier.close()

    
    return
    


def FLAG_DATA_CONCATENATE():

    if REFLAGGING_DATA == True:
        print('START FLAGGING')


    for k in range(len(DATA_band_name)):
        
        for filenames in glob.glob(PROCESS_DIR+'/'+DATA_band_name[k]+'/*.fits'):
                    
            wlmin = DATA_band_min[k]*1E-6
            wlmax = DATA_band_max[k]*1E-6
                    
        
            with fits.open(filenames, memmap=False, mode='update') as fichier:
    
                # print(filenames)
                # name = np.zeros(len(data), dtype=object)
                # for i in range(len(data)):
                #     name[i] = data[i].name
    
                # data=fits.open(filenames, mode='update', memmap=False)
    
                # Read data
                name_HDU = np.array([fichier[t].name for t in range(len(fichier))])
                
                index_vis = np.where(name_HDU=='OI_VIS2')[0]
                
                index_wvl = np.where(name_HDU=='OI_WAVELENGTH')[0]
                
                wavel_shape = np.array([len(fichier[t].data['EFF_WAVE']) for t in index_wvl])

                for i in range(len(index_vis)):
                
                    VIS2     = fichier[index_vis[i]].data['VIS2DATA']
                    VIS2ERR  = fichier[index_vis[i]].data['VIS2ERR']
                    
                    
                    isnan = VIS2ERR != VIS2ERR
                
                    # print(isnan)

                    VIS2ERR[isnan] = 1E99
                    
                                        
                    index = index_wvl[np.where(np.shape(VIS2)[1] == wavel_shape)[0][0]]
                                        
                    WLEN = np.array([fichier[index].data['EFF_WAVE']]*np.shape(VIS2)[0])
                # VIS2FLAG = data['OI_VIS2'].data['FLAG']
            
                    vis    = np.sqrt(np.abs(VIS2)) * np.sign(VIS2)
                    viserr = 0.5* VIS2ERR / (vis+(vis==0));
                    
                    if REFLAGGING_DATA == True:
                        
    
                        
                        flag = (WLEN < wlmax)      &\
                               (WLEN > wlmin)      &\
                                (vis > 0. - viserr) &\
                                (vis < 1. + viserr) &\
                                (vis > -0.1)         &\
                                (vis < 1.1)          &\
                                (viserr > 0)         &\
                                (viserr < 0.1)    #   &\
    
                    else:
                        
                        flag = (WLEN < wlmax)      &\
                               (WLEN > wlmin)

                        
                    flag = np.logical_or(~flag,fichier[index_vis[i]].data['FLAG'])
                    try:
                        fichier[index_vis[i]].data['FLAG'] = flag    
                    except:
                        fichier[index_vis[i]].data['FLAG'] = flag.reshape((len(flag),))
                    
                fichier.flush()
                fichier.close()

    
    return