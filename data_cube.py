# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:35:42 2021

@author: jdrevon
"""
from astropy.io import fits 
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def image_to_data_cube(images, x_axis , y_axis, wavel, PATH, CUBE_NAME):
        
    # This routine returns a FIT file from series of 2D images
        
    # I : M images of size NxK
    # I : x_axis of size N
    # I : y_axis of siez K
    # I : wavel of size M
    
    # 0 : MxNxK image with intensity ratio
    
    size_x = len(x_axis)
    size_y = len(y_axis)
    
    try:
        dimension = len(images)
        cube = np.zeros((dimension,size_x,size_y))
    except:
        dimension = 1
    
    cube = np.zeros((dimension,size_x,size_y))
    

    for i in range(dimension):
        
        try:
            cube[i,:,:] = images[i]
        except:
            cube[i,:,:] = images
    
    primary_hdu = fits.PrimaryHDU(x_axis)
    image_hdu = fits.ImageHDU(images)
    image_hdu2 = fits.ImageHDU(wavel)
    hdu_new = fits.HDUList([primary_hdu,image_hdu])
    hdu_new.append(image_hdu2)
    hdu_new.writeto(PATH+ CUBE_NAME + '.fits',overwrite=True)
    hdu_new.writeto(PATH+ CUBE_NAME +'.fits',overwrite=True)

    return 



