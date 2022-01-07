#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 10:58:06 2021

@author: jdrevon
"""

############ DIRECTORIES:
    
# Path of the folder containing all the OIFITS final calibrated files (LM+N):

#DATA_DIR = "C:/Users/jdrevon/Desktop/Github/RHAPSODY/TEST/NOMEANBCD"

DATA_DIR = "C:/Users/jdrevon/Desktop/R_for/DATA"

# FOLDER USED TO STOCK THE DATA AND DISPLAY THE RESULTS

PROCESS_DIR = "C:/Users/jdrevon/Desktop/Github/RHAPSODY_Rfor/TEST/RFor"

############ PARALLELIZATION:

# Number of telescopes used

ntelescope = 4

############ PARALLELIZATION:
    
nprocs = 8 # Nbr of cores 

############ DATA ERROR BARS:
    
ERROR_SUP_LM = 0.03

ERROR_SUP_N  = 0.03

# Does the error bar on visibilities seems underestimated? 
#The real error bars are replaced by np.sqrt(V2_err**2+ERROR_SUP**2) during the computation! 
# The OIFITS files are not modified

############ DATA WVL ranges

MATISSE_L_BAND_min = 2.96   #µm
MATISSE_L_BAND_max = 4.0   #µm

MATISSE_N_BAND_min = 8.5   #µm
MATISSE_N_BAND_max = 12.5   #µm

########### RINGS:

size_ring_LM = 2   # size_ring_LM [mas] : constant width of a ring (outter diameter-inner diameter) in LM band
size_ring_N  = 5   # size_ring_LM [mas] : constant width of a ring (outter diameter-inner diameter) in N band

NBR_ring_LM  = 56   # NBR_ring_LM [#] : number of rings to put in the model in LM band
NBR_ring_N   = 60   # NBR_ring_N  [#] : number of rings to put in the model in N band

init = 'G' # G for Gaussian like profile, PW for power-law like profile

# For PW profile:

alpha_LM     = .2   # alpha_LM [#] : Power law of the intensity profile for LM band (1/r**(alpha))
alpha_N      = .2   # alpha_N  [#] : Power law of the initial guess intensity profile for N band (1/r**(alpha))

# For G profile:
    
sigma_LM     = 5    # sigma_LM [mas] : Standard deviation of the intensity profile for LM band 
sigma_N      = 10   # sigma_N  [mas] : Standard deviation of the intesnity profile for N band 

########### Method of Regularization: 

# Options : Total Variation (TV) or Quadratic Smoothness (QS)
# Recommanded : Total Variation

#ftot(λ)=χ2(λ)+μfprior(λ)

REG_method = 'TV' #fprior

########### Value of the Hyperparameter: 

HP = [1E0,1E1,1E2,1E3,2E3,4E3,6E3,8E3,1E4,1E5] # µ

########### Fitting Parameters:
    
# WARNING: Only the COBYLA method has been fully tested for this code.
# Please see the different method available here: https://lmfit.github.io/lmfit-py/fitting.html
 
fitting_method  = 'COBYLA' 
tolerance       = 5E-3 #5E-4 fitting and parameter tolerance see lmfit documentation for further explanation (the lower the tolerance value the higher the precision but the longer the computer time)
max_iterations  = 1E6  # Maximum of iterations before stopping for non-convergence


########### Set the resolution of the modeled curve for the plot:

# L-band    

model_q_min_LM   = None # [rad^-1] If None the model boundaries takes the minimum of the spatial frequency covered by the instrument
model_q_max_LM   = None # [rad^-1] If None the model boundaries takes the maximum of the spatial frequency covered by the instrument
model_q_LM       = 1000   # Number of points in the model  

model_rho_min_LM   = None # [rad^-1] If None the intensity profile starts at 0 (the stallar center)
model_rho_max_LM   = None # [rad^-1] If None the intensity profile ends at the model edges
model_rho_LM       = 1000  # Number of points in the model  

# N-band    

model_q_min_N   = None # [rad^-1] If None the model boundaries takes the minimum of the spatial frequency covered by the instrument
model_q_max_N   = None # [rad^-1] If None the model boundaries takes the maximum of the spatial frequency covered by the instrument
model_q_N       = 1000    # Number of points in the model  

model_rho_min_N   = None # [rad^-1] If None the intensity profile starts at 0 (the stallar center)
model_rho_max_N   = None # [rad^-1] If None the intensity profile ends at the model edges
model_rho_N       = 1000  # Number of points in the model  


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

# This second routine aims to gives you the choice to run the fitting process or not. If not, the routine will provide you as a result the initial guess
# parameter entered by the user.

FITTING = True

#3 :  Instrument spectra normalization with a black body. 
# This routine will normalize the instrument spectra by a blackbody spectrum.
# In order to do it the routine needs the black body temperature, the estimated stellar radius in astronomical units and the distance in parsec of the object.

BB_norm = True 

BB_temperature = 2700 # Kelvins
distance_target = 632 #pc
stellar_radii = 3 # AU


