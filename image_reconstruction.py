# -*- coding: utf-8 -*-
"""
Created on Mon May  2 19:09:19 2022

@author: jdrevon
"""

from pretty_table_reading import READ_PRETTY_TABLE
import numpy as np
from rhapsody_init import ERROR_SUP, DATA_DIR, DATA_band_name, PROCESS_DIR
from initialisation_rings import rings
from fits_reading_dico import OIFITS_READING, OIFITS_SORTING
from stock_dico_values import stock_V2_from_dico
from data_copy import COPY_DATA
from model_visibilities import V_ring
from scipy.interpolate import interp1d , interp2d
import matplotlib.pyplot as plt
from numpy import cos as cos
from numpy import sin as sin
import scipy
from imaging_model import imaging_model
from scipy.ndimage import rotate 
from plot_intensity_rec_model import plot_intensity_model_image
from data_cube import image_to_data_cube
def coord_rotation(qu, qv, ang_tmp):
    
    ang             = ang_tmp*np.pi/180
    
    new_qu          = qu*cos(ang)+qv*sin(ang)
    new_qv          = -qu*sin(ang)+qv*cos(ang)

    new_q = np.sqrt((new_qu)**2 + (new_qv)**2)

    return new_qu, new_qv, new_q    


def coord_inclination(qu, qv, inclination_tmp):
    
    inclination     = inclination_tmp*np.pi/180
    
    new_qu          = qu
    new_qv          = qv*cos(inclination)

    new_q = np.sqrt((new_qu)**2 + (new_qv)**2)

    return new_qu, new_qv, new_q    

    
def inv_coord_rotation(qu, qv, ang_tmp):
    
    ang             = ang_tmp*np.pi/180
    
    new_qu          = qu*cos(np.pi/2-ang)+qv*sin(np.pi/2-ang)
    new_qv          = -qu*sin(np.pi/2-ang)+qv*cos(np.pi/2-ang)

    new_q = np.sqrt((new_qu)**2 + (new_qv)**2)

    return new_qu, new_qv, new_q    


def inv_coord_inclination(qu, qv, inclination_tmp):
    
    inclination     = inclination_tmp*np.pi/180
    
    new_qu          = qu
    new_qv          = qv*cos(inclination)

    new_q = np.sqrt((new_qu)**2 + (new_qv)**2)

    return new_qu, new_qv, new_q    



# # COPY THE DATA IN THE RIGHT FOLDER 

# COPY_DATA(DATA_DIR, PROCESS_DIR)

# # SORT DATA in FOLDERS based on V2 and T3 flags:
    
# OIFITS_SORTING()
    
# # READ ALL THE DATA and put them in dicos in order to manipulate the data easily

# OIFITS_TOT = OIFITS_READING()
    

# # We regroup all the usefull information that we want in simple masked array to handle easily the data for the fitting

# wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR, U_DATA, V_DATA = stock_V2_from_dico(OIFITS_TOT)

# # Let's correct the ERROR BARS here:

# for i in range(len(DATA_band_name)):
    
#     V2_DATA_ERR[i] = np.sqrt(V2_DATA_ERR[i]**2+ERROR_SUP[i]**2)
    
# qu_DATA = U_DATA/wavel_DATA
# qv_DATA = V_DATA/wavel_DATA

# # Rings initialization
    
# diam_inner_ring, diam_outter_ring, I_norm, flux, flux_max = rings()


def image_reconstruction(PATH_OUTPUT_FIT_RES, PATH_OUTPUT_IMAGE_REC, band_name, diam_inner_ring, diam_outter_ring, bounds):
    
    PATH_flux = PATH_OUTPUT_FIT_RES+'/fit_flux_ratio_'+band_name+'_band.dat'
    
    header_flux,data_flux = READ_PRETTY_TABLE(PATH_flux, 1)
    
    wavel = data_flux.T[0]
    flux  = data_flux.T[1:].T
    
    PATH_angle = PATH_OUTPUT_FIT_RES+'/angle_'+band_name+'_band.dat'
    header_angle,data_angle = READ_PRETTY_TABLE(PATH_angle, 1)
    
    angle_all =  data_angle.T[1:][0]
    
    PATH_inc = PATH_OUTPUT_FIT_RES+'/inclination_'+band_name+'_band.dat'
    header_inc,data_inc = READ_PRETTY_TABLE(PATH_inc, 1)
    
    inc_all = data_inc.T[1:][0]
    
    
    # Pre-Computing for each bandwidth the associated modeled visibilities for each rings
        
                        
    # qu_interp = np.array([0]+np.logspace(4, 9, 1500).tolist())
    # qv_interp = np.array([0]+np.logspace(4, 9, 1500).tolist())

    qu_interp = np.array([0]+np.logspace(4, 9, 2000).tolist())
    qv_interp = np.array([0]+np.logspace(4, 9, 2000).tolist())

    V_tmp = [V_ring(qu_interp, qv_interp, diam_inner_ring[j],diam_outter_ring[j]) for j in range(len(diam_outter_ring))]            
    q_interp=np.sqrt(qu_interp**2+qv_interp**2)
    V_model_TOT = V_tmp
    
    stock_images = []    
    
    print('START PLOTTING INTENSITY RADIAL PROFILES AND EQUIVALENT INCLINED/ROTATED 2D IMAGES')

    for m in range(len(wavel)):

        print('Wavelength n°%i/%i'%(m+1,len(wavel)))

        res_F = flux[m]

        inc =  inc_all[m]
        ang  =  angle_all[m]
    
    
        V_model= V_model_TOT
        
        CF = np.sum([V_model[i]*res_F[i] for i in range(len(res_F))],axis=0)     
        F_tot = np.sum([res_F[k] for k in range(len(res_F))],axis=0)
        
        Vis_tmp = CF/F_tot
    
        q, image  = imaging_model(q_interp, Vis_tmp, max(q_interp))
        
        qu = q
        qv = q
        
        qu_inc, qv_inc, q_inc = coord_inclination(qu, qv, inc)
        f_inc = scipy.interpolate.interp2d(qu, qv, image)
        image_inc = f_inc(qu_inc,qv_inc)
        
        qu_rot, qv_rot, q_rot = inv_coord_rotation(qu_inc, qv_inc, ang)
        image_rot = rotate(image_inc,ang, reshape= False)

        
        # plt.figure()
        # plt.imshow(np.abs(image_rot), extent=[min(qu),max(qu),min(qv),max(qv)])
        
        from image_to_FFT import FFT_to_image
        
        x, image = FFT_to_image(q, image_rot)
        
        # print(x)
                
        image = np.abs(np.fft.ifftshift(image))/np.amax(np.real(image))
        
        plot_intensity_model_image(x, image, wavel[m], bounds, SAVE_OUTPUT=PATH_OUTPUT_IMAGE_REC)
        
        func = interp2d(x, x, image)
        
        x_new = np.linspace(-bounds,bounds, 512)
        
        cropped = func(x_new, x_new)
        
        # print(np.shape(cropped))
        ########################################         
                       
        # res = cv2.resize(image, dsize=(512,512), interpolation=cv2.INTER_CUBIC)
        # x_new = np.linspace(min(x),max(x),512)
        stock_images.append(cropped)

    # print(np.shape(stock_images))
    image_to_data_cube(stock_images, x_new , x_new, wavel, PATH_OUTPUT_FIT_RES, 'image_datacube_inc_rot')
    print('END PLOTTING INTENSITY RADIAL PROFILES AND EQUIVALENT INCLINED/ROTATED 2D IMAGES')

        
        
        # x_max  = max(x)
    
        # plt.figure()
        # plt.imshow(np.abs(np.fft.ifftshift(image))/np.amax(np.real(image)), extent = [-x_max, x_max, -x_max, x_max])

    return     