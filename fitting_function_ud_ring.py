#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 15:01:53 2021

@author: jdrevon
"""

import numpy as np

from rhapsody_init import REG_method, max_iterations, tolerance, nprocs, DATA_band_name
from lmfit import minimize as minim
from lmfit import Parameters,fit_report
from model_visibilities import V_uniform, V_ring
from check_intensity import check_intensity
from plot_V2_obs_model import plot_V2_obs_model
from f_prior import total_variation, quad_smoothness
from plot_hist import plot_histogram
from prettytable import PrettyTable
import multiprocessing as mp

def short_model(params,q_UD_fitting, V2_DATA, V2_DATA_ERR, HP, V_model, REG_method, intensity_init, flux_init):
    
    
    nb_ring = len(flux_init)
    res_diam_outer_ring         = [params['diam_outer_ring'+str(i)].value for i in range(nb_ring)]
    res_diam_inner_ring         = [params['diam_inner_ring'+str(i)].value for i in range(nb_ring)]
    res_F_tmp                   = [params['flux_ftot'+str(i)].value for i in range(nb_ring)]

    res_F, I_tot_norm = check_intensity(res_F_tmp, res_diam_outer_ring, res_diam_inner_ring, intensity_init, flux_init)

    CF = np.sum([V_model[i]*res_F[i] for i in range(nb_ring)],axis=0)

     
    F_tot = np.sum([res_F[k] for k in range(nb_ring)],axis=0)
    Vis2 = (CF/F_tot)**2
   
    # plt.figure()
    # plt.scatter(q_UD_fitting,Vis2,s=2)
    
    if REG_method == 'TV':

        f_prior = total_variation(HP, I_tot_norm, np.array(res_diam_outer_ring)/2)    

    elif REG_method == 'QS':

        f_prior = quad_smoothness(HP, I_tot_norm)    
        
    residuals = (Vis2-V2_DATA)/V2_DATA_ERR 
    chi2_tot   = np.sum(residuals**2) + f_prior

    return chi2_tot

def fitting_function(i, wavel_ALL, wavel_DATA, HP, diam_outter_ring, \
                     params2, \
                         q_UD_HR, q_DATA, \
                             intensity_profile, \
                                 flux_ratio, \
                                     V2_DATA, V2_DATA_ERR, V_model, \
                                         PATH_OUTPUT_VIS, PATH_OUTPUT_HIST):
    
    nb_ring = len(flux_ratio)    
    cond= wavel_DATA == wavel_ALL[i]
    q_UD_fitting      = np.array(q_DATA[cond])  
    V2_DATA     = np.array(V2_DATA[cond])
    V2_DATA_ERR = np.array(V2_DATA_ERR[cond])                         


    res2 = minim(short_model, params2, args=(q_UD_fitting, V2_DATA, V2_DATA_ERR, HP, V_model[i], REG_method, intensity_profile, flux_ratio), method='COBYLA',max_nfev=max_iterations, tol = tolerance) # 5E-4 OK

    print('WORK ON %.3f OVER'%(wavel_ALL[i]*1E6))
    
    res_diam_outer_ring = [res2.params['diam_outer_ring'+str(k)].value for k in range(nb_ring)]
    res_diam_inner_ring = [res2.params['diam_inner_ring'+str(k)].value for k in range(nb_ring)]
    res_F_tmp = [res2.params['flux_ftot'+str(k)].value for k in range(nb_ring)]
    
    # This function has been already included in the computation of the chi2 but have to be also execute in the output to take into consideration the different changes made on the flux values.
    # Here we say that wa cannot have an intensity higher than the central uniform disk. Which is in the continuity of our initial condition.
    
    
    res_F, I_tot_norm = check_intensity(res_F_tmp, res_diam_outer_ring, res_diam_inner_ring, intensity_profile, flux_ratio)


    # Indeed, it is an easy task to compute I_norm (the normalized intensity profil with respect to the central uniform disk value) from the flux since we control the geometry of our model.
    # However in order to gain some computation time and since we have already computed I_norm, we set both as inputs instead of only one of the two.

    # Then based on this new flux ratio for each rings, we compute the coherent flux in order to compute the visibilities

    CF_fitting = np.sum([V_model[i][k]*res_F[k]\
                for k in range(nb_ring)], axis=0)

            
    F_tot = np.sum([res_F[k] for k in range(nb_ring)],axis=0)
    V2_fitting = (CF_fitting/F_tot)**2

    # Here based on the final fitting parameters for this model at this wl, for display purpose, we want to increase the resolution of the model plots
    # Then we compute a last time the visibility curves.

    CF = np.sum([V_ring(q_UD_HR,res_diam_inner_ring[k],res_diam_outer_ring[k])*res_F[k]\
                for k in range(nb_ring)],axis=0)

    V2_model_HR = (CF/F_tot)**2


    # PLOT of the V2 with higher resolution for the model displayed in linear scale            

    plot_V2_obs_model(q_UD_fitting,V2_DATA, q_UD_HR, V2_model_HR, wavel_ALL[i]*1E6, PLOT = False,\
                      SAVE_OUTPUT = PATH_OUTPUT_VIS, logy = True, y_obs_err = V2_DATA_ERR,\
                          xlim_min = 1E6, xlim_max = 5E7, ylim_min = 1E-7, ylim_max = 1E0)           

    # PLOT of the V2 with higher resolution for the model displayed in logarithmic scale            

    plot_V2_obs_model(q_UD_fitting,V2_DATA, q_UD_HR, V2_model_HR, wavel_ALL[i]*1E6, PLOT = False,\
                      SAVE_OUTPUT = PATH_OUTPUT_VIS, y_obs_err = V2_DATA_ERR,\
                          xlim_min = 0, xlim_max = 5E7, ylim_min = -0.1, ylim_max = 1.05)           

   
    # Here we set an histogram of the residuals weighted by the error in order to check if we have a Gaussian like distribution to quantify the goodness of fit in addition to the chi2 value
    # However knowing that we use a Bayesian approach the more we increase the prior impact on the fit, the more we will observe an expected deviation from the Gaussian like distribution
   
    plot_histogram(V2_DATA,V2_DATA_ERR,V2_fitting, wavel_ALL[i]*1E6, PLOT = False, SAVE_OUTPUT = PATH_OUTPUT_HIST)
    
    
    # Now we manually compute the chi2 value and chi2_reduced of the best model fitting at this wavelength taking into account the regularization 
    
    # CHI2:

    N = len(V2_DATA)
    chi2_tmp = np.sum(((V2_DATA-V2_fitting)/V2_DATA_ERR)**2)
    # chi2[i]=chi2_tmp

    # CHI2_REDUCED:

    if nb_ring>1 :
        D_number = N-nb_ring 
        chi2_red= np.abs(chi2_tmp/D_number)
        
    else:
        None
        
    # Regularization:
        
    if REG_method == 'TV':

        f_prior_tmp = total_variation(HP, I_tot_norm, np.array(res_diam_outer_ring)/2)    

    elif REG_method == 'QS':

        f_prior_tmp = quad_smoothness(HP, I_tot_norm)    
           
    # CHI2_TOT
    
    chi2_tot_tmp = chi2_tmp + f_prior_tmp

    # INTENSITY PROFILE
    
    res_F_mod = np.array(res_F)/np.sum(res_F) 
    res_F = np.round(res_F_mod,6)
    I_tot = I_tot_norm

    # print(fit_report(res2))

    # print('WORK ON %.3f OVER'%wavel_ALL[i])

    
    return  wavel_ALL[i], res2, chi2_tmp, chi2_red, f_prior_tmp, chi2_tot_tmp, res_F, I_tot



def UD_modeling(wavel_DATA, q_DATA, V2_DATA, V2_DATA_ERR,\
                         q_UD_HR, diam_outter_ring, diam_inner_ring, flux_ratio, I_norm, flux_max,\
                                     HP, V_model,\
                                         PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES, PLOT):
         

        

    # INITIALIZATION OF THE MINIMIZATION ROUTINE PARAMETER'S 

    name_var_ring = [['diam_outer_ring'+str(i) for i in range(len(diam_outter_ring[k]))] for k in range(len(DATA_band_name))]
    name_var_ring_bef = [['diam_inner_ring'+str(i) for i in range(len(diam_inner_ring[k]))] for k in range(len(DATA_band_name))]
    name_var_F = [['flux_ftot'+str(i) for i in range(len(flux_ratio[k]))] for k in range(len(DATA_band_name))]

    variation_flux = [np.ones(len(diam_outter_ring[k])).astype('bool') for k in range(len(DATA_band_name))]
    for k in range(len(DATA_band_name)): variation_flux[k][0] = False

    params_TOT = []
    for k in range(len(DATA_band_name)):
        
        params2 = Parameters()        
        [params2.add(name_var_ring[k][j],value = diam_outter_ring[k][j], vary= False) for j in range(len(diam_outter_ring[k]))]
        [params2.add(name_var_ring_bef[k][j],value = diam_inner_ring[k][j], vary= False) for j in range(len(diam_inner_ring[k]))]
        [params2.add(name_var_F[k][j],value = flux_ratio[k][j], min=0, max=flux_max[k][j], vary= variation_flux[k][j]) for j in range(len(diam_outter_ring[k]))]

        # Indeed the more we get further from the star, the more the rings have low probability to have higher intensity than the central uniform disk using this maximum boundary
        # This will also hardly depend on the initial intensity profile that we use. At the moment only 1/r^alpha are available thus this is impossible 
        # to have high intensity rings far from the star.
        
        params_TOT.append(params2)

    
    # OK LET'S GO FIT: DON'T STOP ME NOW !
    
    for k in range(len(DATA_band_name)):

        list_wavel = np.unique(wavel_DATA[k])
        
        print('STARTING FITTING on %s band'%DATA_band_name[k])
    
        pool = mp.Pool(processes=nprocs)    
        result_parallel = pool.starmap(fitting_function, [(i, list_wavel, wavel_DATA[k], HP, diam_outter_ring[k], \
                             params_TOT[k], q_UD_HR[k], q_DATA[k], \
                                     I_norm[k], flux_ratio[k], \
                                             V2_DATA[k], V2_DATA_ERR[k], V_model[k], \
                                                 PATH_OUTPUT_VIS[k], PATH_OUTPUT_HIST[k]) for i in range(len(list_wavel))])
    
        pool.close()
        pool.join()
    
    
    
        list_wavel, resultat, chi2, chi2_red, f_prior, chi2_tot, res_F, I_tot = np.array(result_parallel).T

        print('END FITTING on %s band'%DATA_band_name[k])
        
        # TABLE OF RESULTS
    
    
        # INITIALIZATION OF THE RESULTS TABLE OF THE FITTING FLUX FOR EACH COMPONENTS OF THE MODEL: DISPLAY ==> FLUX vs RADIUS  
        
        list_wavel = list_wavel*1E6
        
        table_res_fit = PrettyTable()
        table_res_fit.field_names = ["WAVEL"] + [str((np.array(diam_outter_ring[k])[i]/2)) for i in range(len(diam_outter_ring[k]))]
    

    
        for ind in range(len(list_wavel)):        
            table_res_fit.add_row([list_wavel[ind]] + np.round(res_F[ind],6).tolist())
    
        # TABLE OF FLUX
    
        print('GENERATING FLUX TABLE for the %s band'%DATA_band_name[k])
    
        with open(PATH_OUTPUT_FIT_RES[k]+'fit_flux_ratio_%s_band.dat'%DATA_band_name[k], 'w') as f: f.write(str(table_res_fit))
    
        # TABLE OF INTENSITY
        
        print('GENERATING INTENSITY TABLE for the %s band'%DATA_band_name[k])
        
        data_intensity = PrettyTable()
        column_names_I = ["rho_ext"]+[str(np.round(list_wavel[i],5)) for i in range(len(I_tot))]
        data_intensity.add_column(column_names_I[0], (np.array(diam_outter_ring[k])/2).tolist())
        for i in range(len(I_tot)):
            data_intensity.add_column(column_names_I[i+1],I_tot[i])
        with open(PATH_OUTPUT_FIT_RES[k]+'intensity_%s_band.dat'%DATA_band_name[k], 'w') as f: f.write(str(data_intensity))
        
        # ALL LMFIT REPORT
    
        print('GENERATING FITTING REPORT for the %s band'%DATA_band_name[k])
    
    
        with open(PATH_OUTPUT_FIT_RES[k]+ 'resultat_fitting_rings_%s_band.txt'%DATA_band_name[k], "w+") as att_file:
            for ind in range(len(list_wavel)):
                att_file.write(str(list_wavel[ind])+ ' µm' + '\n')
                att_file.write(fit_report(resultat[ind]) + '\n')
    
        # ALL CHI2 REPORT
    
        print('GENERATING CHI2 TABLE for the %s band'%DATA_band_name[k])
    
            
        data_chi2 = PrettyTable()
        data_chi2.field_names = ['Wavel', 'chi2', 'chi2_red', 'f_prior_'+ REG_method, 'chi2_tot']
        for j in range(len(chi2)):
            data_chi2.add_row([list_wavel[j]]+ [chi2[j]] + [chi2_red[j]] + [f_prior[j]] + [chi2_tot[j]])
        with open(PATH_OUTPUT_FIT_RES[k]+'chi2_%s_band.dat'%DATA_band_name[k], 'w') as f: f.write(str(data_chi2))
    
    return


### Optimisation


    # import fileinput
    # import cProfile
    
    # pr = cProfile.Profile()
    # pr.enable()

    # pr.disable()
    # filename = 'profile.prof'  # You can change this if needed
    # pr.dump_stats(filename)
    # import sys
    # from pstats import Stats
    
    # my_stat = Stats('profile.prof', stream=sys.stdout)
    # my_stat.sort_stats('time').print_stats(10)









