# -*- coding: utf-8 -*-
"""
Created on Mon May  2 22:30:41 2022

@author: jdrevon
"""


from A1_mas_to_rad import mas_to_rad, rad_to_mas
import numpy as np
from imaging_model import imaging_model
import matplotlib.pyplot as plt
import numpy.fft as fft
from numpy import cos as cos
from numpy import sin as sin
import scipy

def coord_change(qu, qv, inclination_tmp, ang_tmp):
    
    ang             = ang_tmp*np.pi/180
    inclination     = inclination_tmp*np.pi/180
    
    new_qu          = (qu*cos(ang)+qv*sin(ang))*cos(inclination)
    new_qv          = -qu*sin(ang)+qv*cos(ang)

    new_q = np.sqrt((new_qu)**2 + (new_qv)**2)

    return new_qu, new_qv, new_q    
    

def image_to_FFT(x, image):
    
    #the image should be squared and centered
    
    #### COMPUTE THE FFT

    nb_points = len(x)
    q_delta = 1/(2*mas_to_rad(max(x)))
    q_max   = nb_points*q_delta/2
    
    q = np.linspace(-q_max,q_max,nb_points) # rad^-1
    
     
    FFT = fft.fftshift(fft.fft2(image))

    # plt.figure()
    # plt.imshow(np.abs(FFT), extent=[-q_max,q_max,-q_max,q_max])

    return q, FFT


def FFT_to_image(q, FFT):
    
    q_delta = max(q)*2/len(q)
    champ  = 1/q_delta
    champ_mas = rad_to_mas(champ)

    x = np.linspace(-champ_mas/2,champ_mas/2, len(q))    
    
    image = fft.ifft2(FFT)
    
    return x, image


# champ       = 80 # mas 
# diam         = 10 # mas
# nb_points    = 100 
# rho          = np.linspace(0,champ/2,nb_points)

# I_profile              = np.zeros(nb_points)
# I_profile[rho<=diam/2] = 1


# xv, image = imaging_model(rho, I_profile, champ/2)

# ##############################


# x_mod, y_mod, x_tot_mod = coord_change(xv, xv, 0, 0)

# f_image = scipy.interpolate.interp2d(x_mod, y_mod, image)

# image_mod = f_image(xv,xv)


# # plt.figure()
# # plt.imshow(image_mod,extent=[min(xv),max(xv),min(xv),max(xv)])

# image_rotated = function(xv_mod,xv_mod)


# I_function = scipy.interpolate.interp1d(rho, I_profile)



# rho_x = np.linspace(-max(rho),max(rho),len(rho),endpoint=False)
# rho_y = np.linspace(-max(rho),max(rho),len(rho),endpoint=False)

# x_transform, y_transform, x_tot = coord_change(rho_x, rho_y, 75, 20)

# x_transform_mesh,y_transform_mesh = np.meshgrid(x_transform[25:-25], y_transform[25:-25])

# I_new = I_function(np.sqrt(x_transform_mesh**2+y_transform_mesh**2))


# xv_new_mesh, xv_new_mesh, x_tot = coord_change(xv_mesh, xv_mesh, 45, 20)


# plt.figure()
# plt.imshow(image**(1/2), extent=[min(xv),max(xv),min(xv),max(xv)])
# plt.title('Image')

# interp_coord_image = scipy.interpolate.griddata((xv_mesh.ravel(),xv_mesh.ravel()), image.ravel(), (xv_new_mesh.ravel(),xv_new_mesh.ravel()))


# plt.figure()
# plt.imshow(interp_coord_image**(1/2), extent=[min(xv),max(xv),min(xv),max(xv)])
# plt.title('Image after mod.')

# plt.figure()
# plt.imshow(image_rotated**(1/2), extent=[min(xv),max(xv),min(xv),max(xv)])
# plt.title('Image modifiée')


# q, FFT = image_to_FFT(xv, image)


# qu_tmp = q
# qv_tmp = q

# qu_new_tmp,qv_new_tmp,q_new_tmp= coord_change(qu_tmp, qv_tmp, 45, 45)


# qu,qv = np.meshgrid(qu_tmp, qv_tmp)

# qu_new, qv_new = np.meshgrid(qu_new_tmp, qv_new_tmp)


# fourier_mod = scipy.interpolate.griddata((qu.ravel(),qv.ravel()), FFT.ravel(), (qu_new,qv_new))


# plt.figure()
# plt.imshow(np.abs(FFT)**(1/2), extent=[min(q),max(q),min(q),max(q)])
# plt.title('Before change of coordinates')


# plt.figure()
# plt.imshow(np.abs(fourier_mod)**(1/2), extent=[min(q),max(q),min(q),max(q)])
# plt.title('After change of coordinates')


# x_new, image_new = FFT_to_image(q, fourier_mod)

# plt.figure()
# plt.imshow(np.abs((image_new))**(1/2), extent=[min(x_new),max(x_new),min(x_new),max(x_new)])
# plt.title('Image after change of coordinates')
