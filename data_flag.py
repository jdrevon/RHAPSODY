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
                WLEN = np.array([fichier['OI_WAVELENGTH'].data['EFF_WAVE']]*np.shape(VIS2)[0])
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

                        
                flag = ~flag
                fichier['OI_VIS2'].data['FLAG'] = flag    
                
                # print(data['OI_VIS2'].data['FLAG'])
                
                CPERR  = fichier['OI_T3'].data['T3PHIERR']
                WLEN = np.array([fichier['OI_WAVELENGTH'].data['EFF_WAVE']]*np.shape(CPERR)[0])
            
                flag3 = (WLEN < wlmax) &\
                       (WLEN > wlmin) &\
                       (CPERR > 0.)     &\
                       (CPERR < 20.)
                flag3 = ~ flag3
            
                fichier['OI_T3'].data['FLAG'] = flag3
                
                fichier.flush()
                fichier.close()

    
    return
    
