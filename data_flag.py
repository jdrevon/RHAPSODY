# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:52:31 2020

@author: jdrevon
"""

# Import stuff
import astropy
from astropy.io import fits
import numpy as np
from shutil import copyfile
from os import walk
import glob
from rhapsody_init import MATISSE_DIR, MATISSE_L_BAND_min, MATISSE_L_BAND_max, MATISSE_N_BAND_min, MATISSE_N_BAND_max

def FLAG_DATA():

    for filenames in glob.glob(MATISSE_DIR+'*.fits'):
    
        # print(filenames)
    
        if 'IR-N' in filenames:
            
            wlmin = MATISSE_L_BAND_min*1E-6
            wlmax = MATISSE_L_BAND_max*1E-6
        
        elif 'IR-LM' in filenames:
            
            wlmin = MATISSE_N_BAND_min*1E-6
            wlmax = MATISSE_N_BAND_max*1E-6
    
        with fits.open(filenames, memmap=False, mode='update') as data:

        # data=fits.open(filenames, mode='update', memmap=False)
        
            WLEN = data['OI_WAVELENGTH'].data['EFF_WAVE']
        
            # Read data
            
            VIS2     = data['OI_VIS2'].data['VIS2DATA']
            VIS2ERR  = data['OI_VIS2'].data['VIS2ERR']
            # VIS2FLAG = data['OI_VIS2'].data['FLAG']
        
            vis    = np.sqrt(np.abs(VIS2)) * np.sign(VIS2)
            viserr = 0.5* VIS2ERR / (vis+(vis==0));
            
            #flag = ~ VIS2FLAG
            flag = (WLEN < wlmax)      &\
                   (WLEN > wlmin)      &\
                   (vis > 0. - viserr) &\
                   (vis < 1. + viserr) &\
                   (vis > -0.1)         &\
                   (vis < 1.1)          &\
                   (viserr > 0)         &\
                   (viserr < 0.1)    #   &\
                  # (VIS2 / VIS2ERR > 3)
        
            #print(flag)
            flag = ~flag
        
            data['OI_VIS2'].data['FLAG'] = flag
            
                        
            # CP     = data['OI_T3'].data['T3PHI']
            CPERR  = data['OI_T3'].data['T3PHIERR']
            # CPFLAG = data['OI_T3'].data['FLAG']
        
            flag3 = (WLEN < wlmax) &\
                   (WLEN > wlmin) &\
                   (CPERR > 0.)     &\
                   (CPERR < 20.)
            flag3 = ~ flag3
        
            data['OI_T3'].data['FLAG'] = flag3
            
            data.flush()
            data.close()

    
    return
    
