#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 10:58:06 2021

@author: jdrevon
"""

############ DIRECTORIES:
    
# Path of the folder containing all the OIFITS final calibrated files:
#Ex: 
    #DATA_DIR = ["C:/Users/jdrevon/Desktop/Article R Scl/R_Scl2/LM",
#         "C:/Users/jdrevon/Desktop/Article R Scl/R_Scl2/N"] #Path for all the different bands


# DATA_DIR = ["C:/Users/jdrevon/Desktop/Presentations/Conference_2022/SPIE/beauty_contest"] #Path for all the different bands

DATA_DIR = ["C:/Users/jdrevon/Desktop/Presentations/Conference_2022/SPIE/MODEL_proceding/Inclined"] #Path for all the different bands


# FOLDER USED TO STOCK THE DATA AND DISPLAY THE RESULTS

PROCESS_DIR = "C:/Users/jdrevon/Desktop/Presentations/Conference_2022/SPIE/MODEL_proceding/Inclined"

# PROCESS_DIR = "C:/Users/jdrevon/Desktop/Presentations/Conference_2022/SPIE/beauty_contest"

# PROCESS_DIR = "C:/Users/jdrevon/Desktop/Presentations/Conference_2022/SPIE/DATA"

############ PARALLELIZATION:

# Number of telescopes used

ntelescope = 4 #4

############ PARALLELIZATION:
    
nprocs = 8 # Nbr of cores 

############ DATA WVL ranges
#Ex:
# DATA_band_name = ['LM','N'] 
# DATA_band_min  = [2.96,8.50]  #µm  
# DATA_band_max  = [4.00,12.50] #µm

DATA_band_name = ['LM'] 
DATA_band_min  = [1]  #µm 
DATA_band_max  = [5] #µm

# DATA_band_name = ['LM', 'N'] 
# DATA_band_min  = [2.8, 8.50]  #µm 
# DATA_band_max  = [4.00, 12.50] #µm


############ DATA:

# Pass the OIFITS data in the REFLAGGING routine to make a more robust flag
# This are the flagin criteria that we are using:
# We do not recommand this on data with very high error bars.

# vis    = np.sqrt(np.abs(VIS2)) * np.sign(VIS2)
# viserr = 0.5* VIS2ERR / (vis+(vis==0));

# flag = (vis > 0. - viserr) &\
#        (vis < 1. + viserr) &\
#        (vis > -0.1)         &\
#        (vis < 1.1)          &\
#        (viserr > 0)         &\
#        (viserr < 0.1)    #   &\

REFLAGGING_DATA = False


# Does the error bar on visibilities seems underestimated? 
#The real error bars are replaced by np.sqrt(V2_err**2+ERROR_SUP**2) during the computation! 
# The OIFITS files are not modified
    
# ERROR_SUP = [0.03,0.03]

ERROR_SUP = [0.02]

# FLUX
# Is the OIFITS_FLUX table is provided in the data? If Yes, please set True, if not set False

OIFITS_FLUX = False


########### RINGS:

# size_ring = [2,5]      # size_ring [mas] : constant angular diameter of a ring in the i-th band.
# NBR_ring  = [56,60] # NBR_ring [#] : number of rings to put in the model for the different bands


size_ring = [2]      # size_ring [mas] : constant angular diameter of a ring in the i-th band.
NBR_ring  = [100]     # NBR_ring [#] : number of rings to put in the model for the different bands

init = 'G' # G for Gaussian like profile, PW for power-law like profile, D for Uniform Disk

# For PW profile:

# alpha  = [.2,.2] # alpha [#] : Power law of the intensity profile for each bands (1/r**(alpha))

alpha  = [.2] # alpha [#] : Power law of the intensity profile for each bands (1/r**(alpha))

# For G profile:
    
# sigma = [3,10]    # sigma [mas] : Standard deviation of the intensity profile for each bands 

sigma = [10]    # sigma [mas] : Standard deviation of the intensity profile for each bands FWHM ~ 2.3*sigma 

# For Uniform Disk 

diam_disk = [10] # [mas] : Diameter of the Uniform Disk

########### Method of Regularization: 

# Options : Total Variation (TV) or Quadratic Smoothness (QS)
# Recommanded : Total Variation

#ftot(λ)=χ2(λ)+μfprior(λ)

REG_method = 'TV' #fprior

########### Value of the Hyperparameter: 

HP = [0, 1E0, 1E1, 1E2, 1E3, 1E4, 1E5] # µ

########### Fitting Parameters:
    
# WARNING: Only the COBYLA method has been fully tested for this code.
# Please see the different method available here: https://lmfit.github.io/lmfit-py/fitting.html
 
fitting_method  = 'COBYLA' 
tolerance       = 5E-4 #5E-4 fitting and parameter tolerance see lmfit documentation for further explanation (the lower the tolerance value the higher the precision but the longer the computer time)
max_iterations  = 1E6  # Maximum of iterations before stopping for non-convergence

# FLATNESS

# The flatness is applied as a change of coordinate

    # new_qu = (qu_DATA*np.cos(ang)-qv_DATA*np.sin(ang))* flatness
    # new_qv = (qu_DATA*np.sin(ang)+qv_DATA*np.cos(ang)) 

    # new_q = np.sqrt((new_qu)**2 +
    #              (new_qv)**2);

# The inclination made on the fourier transform is computed as to be an inclination in the x axis on the image plane.
# The given inclination is then the inclination that the user want in the image plane along the x-axis.

# inc_flag       = True # If you want to activate the inclination option set "True" otherwise set "False" and put 0 as initial guess.
# inc_guess      = [40] # 0 ==> centro-symmetric (no inclination), >0 ==> inclination .
# angle_guess    = [80] # orientation angle of the structure # 0° = no orientation


inc_flag       = True # If you want to activate the inclination option set "True" otherwise set "False" and put 0 as initial guess.
inc_guess      = [50] # 0 ==> centro-symmetric (no inclination), >0 ==> inclination .
angle_guess    = [50] # orientation angle of the structure # 0° = no orientation

########### Set the resolution of the modeled curve for the plot:

model_visibilities = [2000] #Number of points to interpolate the visibilities (the higher the flatness value is, the higher the number of points is required) 

# model_rho_min = [None, None] # [rad^-1] If None the intensity profile starts at 0 (the stallar center)
# model_rho_max = [None, None] # [rad^-1] If None the intensity profile ends at the model edges
# model_rho     = [1000, 1000] # Number of points in the model for each bands


model_rho_min = [None] # [rad^-1] If None the intensity profile starts at 0 (the stallar center)
model_rho_max = [None] # [rad^-1] If None the intensity profile ends at the model edges
model_rho     = [1000] # Number of points in the model for each bands

image_rec_windows = [60] # The maximum radial distance that the user want to display on the image reconstruction and the intensity profile plots

########### OPTIONNAL FEATURES PARAMETERS:
    
# The followings features are not essentials in the RHAPSODY code, however we decided to keep them for users.
# They are still present in the code we kept them, since the allocated computing time for them is negligeable
# Here are the few parameters of those features that you can adjust at your convenience

#1 : Position Angle groups

# This first routine is able to flag the visibility data points and its associated baseline with respect to the Position Angle.
# The user can tell in which different he wants for the in how many equal parts he wants to break down the UV plan in order to define groups.
# For example setting the NBR_GRPS_PA = 2 will devide the |PA| between 0 and 180 in two groups: 0-90 and 90-180.

NBR_GRPS_PA = 4

#2 : To Fit or not to fit... that is the question

# This second routine aims to gives you the choice to run the fitting process or not. If not, the routine will only read the previous data fitted and located in the folders

READING_ONLY = False #READING ONLY

#3 #This routine aims to gives you the choice to run the fitting process with fitting or not. If not, the routine will only run the code with the initial parameters.

FITTING = True

#4 :  Instrument spectra normalization with a black body. 
# This routine will normalize the instrument spectra by a blackbody spectrum.
# In order to do it the routine needs the black body temperature, the estimated stellar radius in astronomical units and the distance in parsec of the object.

BB_norm = True 

BB_temperature = 2700 # Kelvins
distance_target = 632 #pc
stellar_radii = 3 # AU


