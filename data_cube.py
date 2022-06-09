#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:35:42 2021

@author: jdrevon
"""
from astropy.io import fits 
import numpy as np

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
    
    # primary_hdu = fits.PrimaryHDU(x_axis)
    # image_hdu = fits.ImageHDU(images)
    # image_hdu2 = fits.ImageHDU(wavel)

    hdr = fits.Header()
    hdr['OBSERVER']  = 'Julien Drevon'

    hdr['CRPIX1'] = size_x/2
    hdr['CRVAL1'] = 0.0000
    hdr['CDELT1'] = -(min(x_axis)/(size_x/2))#*1E-3*np.pi/(180*3600) 
    hdr['CUNIT1'] = 'mas'

    hdr['CRPIX2'] = size_y/2
    hdr['CRVAL2'] = 0.0000
    hdr['CDELT2'] = -(min(y_axis)/(size_y/2))#*1E-3*np.pi/(180*3600)
    hdr['CUNIT2'] = 'mas'

    hdr['CRPIX3'] = 1
    hdr['CRVAL3'] = min(wavel)
    try:
        hdr['CDELT3'] = np.diff(wavel)[0]
    except:
        hdr['CDELT3'] = 0 # 1 WVL only
    hdr['CUNIT3'] = 'MICRON'
    
    primary_hdu = fits.PrimaryHDU(header=hdr)
    image_hdu = fits.ImageHDU(header=hdr, data=cube[:,:,::-1])
    hdul = fits.HDUList([primary_hdu,image_hdu])
    hdul.writeto(PATH+ CUBE_NAME +'.fits',overwrite=True)

    return 




# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Fri Jul 24 15:25:14 2020

# @author: jdrevon
# """

# from astropy.io import fits 
# import numpy as np
# from scipy.interpolate import interp1d
# import matplotlib.pyplot as plt
# from initial_parameters import DATA_CENTER

# def DUSTY_to_2D(x_intensity,y_intensity,wavelengths,size,boundary,padding):
    
#     # This routine returns a 2D image using the 1D intensity profile of DUSTY:
        
#     # I : x_intensity in arcsec
#     # I : y_intensity in Jy.arcsec^-2 for M wavelengths
#     # I : image size in px (In power of 2)
#     # I : boundary of the image in mas
#     # I : Padding True or False (Go to the next power of 2)
    
#     # 0 : MxNxN image with intensity values in Jy.arsec^-2
#     # 0 : Physical axis values of the squared image in arcsec 
    
#     x_intensity = x_intensity*1E3
    
#     try:
#         dimension = len(y_intensity[0,:])
#         cube = np.zeros((dimension,size,size))
#     except:
#         dimension = 1
#         cube = np.zeros((dimension,size,size,))
    
#     if boundary > max(x_intensity):
#         bounds = max(x_intensity)
#     else :
#         bounds = boundary

#     if padding=='True':
    
#         cube = np.zeros((size*2,size*2,dimension))

#     else : 
#         cube = np.zeros((dimension,size,size))

#     xv = np.linspace(-bounds,bounds,size,endpoint=False)

#     X, Y = np.meshgrid(xv, xv)
#     current_radius = np.sqrt(X**2 + Y**2)
#     cond=np.logical_and(current_radius<=max(x_intensity),current_radius>=min(x_intensity))
    
#     if dimension == 1:

#         y_comp = y_intensity
#         interpol_index = interp1d(x_intensity, y_comp)            
#         #I_min = interpol_index(max(current_radius[cond]))    
#         profilegrid2 = np.zeros(X.shape, float)#*I_min
#         profilegrid2[cond] = interpol_index(current_radius[cond])

#         if padding == 'True':    
#             profilegrid2=np.pad(profilegrid2, ((int(size/2),int(size/2)),(int(size/2),int(size/2))), 'constant') # On padding en haut, en bas, à droite et à gauche
#         cube[0,:,:] = profilegrid2
  

#     else :
        
#         for i in range(dimension):
#             y_comp = y_intensity[:,i]
#             interpol_index = interp1d(x_intensity, y_comp)            
#             #I_min = interpol_index(max(current_radius[cond]))    
#             profilegrid2 = np.zeros(X.shape, float)#*I_min
#             profilegrid2[cond] = interpol_index(current_radius[cond])
        
#             if padding == 'True':    
#                 profilegrid2=np.pad(profilegrid2, ((int(size/2),int(size/2)),(int(size/2),int(size/2))), 'constant') # On padding en haut, en bas, à droite et à gauche
#             cube[i,:,:] = profilegrid2
    
        
#     hdr = fits.Header()
#     hdr['OBSERVER']  = 'Julien Drevon'

#     hdr['CRPIX1'] = 256
#     hdr['CRVAL1'] = 0.0000
#     hdr['CDELT1'] = (-boundary/(size/2))*1E-3*np.pi/(180*3600) 
#     hdr['CUNIT1'] = 'RAD'

#     hdr['CRPIX2'] = 256
#     hdr['CRVAL2'] = 0.0000
#     hdr['CDELT2'] = (boundary/(size/2))*1E-3*np.pi/(180*3600)
#     hdr['CUNIT2'] = 'RAD'

#     hdr['CRPIX3'] = 1
#     hdr['CRVAL3'] = min(wavelengths)
#     hdr['CDELT3'] = np.diff(wavelengths)[0]
#     hdr['CUNIT3'] = 'MICRON'
    
#     primary_hdu = fits.PrimaryHDU(header=hdr)
#     image_hdu = fits.ImageHDU(header=hdr, data=cube)
#     hdul = fits.HDUList([primary_hdu,image_hdu])
#     hdul.writeto(DATA_CENTER+'cube.fits',overwrite=True)
    
#     xv = xv*1E-3
#     if padding == 'True':    
        
#         Ni = size
#         Deltax = 2*bounds
#         deltax = Deltax/Ni  
        
#         xv = np.linspace(-bounds-size/2*deltax,bounds+size/2*deltax,size*2,endpoint=False)*1E-3

#     print('R3_OK')
#     return xv,cube




