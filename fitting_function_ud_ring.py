# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 15:01:53 2021

@author: jdrevon
"""

import numpy as np

from rhapsody_init import model_rho_LM, model_rho_N, REG_method, max_iterations, tolerance, nprocs
from lmfit import minimize as minim
from lmfit import Parameters,fit_report
from model_visibilities import V_uniform, V_ring
from check_intensity import check_intensity
from plot_V2_obs_model import plot_V2_obs_model
from f_prior import total_variation, quad_smoothness
from plot_hist import plot_histogram
from prettytable import PrettyTable
import multiprocessing as mp

import matplotlib.pyplot as plt

def short_model(params,q_UD_fitting, V2_MATISSE, V2_MATISSE_ERR, HP, V_model, REG_method, intensity_init, flux_init):
    
    
    nb_UD = int(len(params)/3)
    res_UD            = [params['diam_UD'+str(i)].value for i in range(nb_UD)]
    res_UD_bef        = [params['diam_UD_bef'+str(i)].value for i in range(nb_UD)]
    res_F_tmp         = [params['flux_ftot'+str(i)].value for i in range(nb_UD)]

    res_F, I_tot_norm = check_intensity(res_F_tmp, res_UD, res_UD_bef, intensity_init, flux_init)

    CF = np.sum([V_model[i]*res_F[i] for i in range(nb_UD)],axis=0)

     
    F_tot = np.sum([res_F[k] for k in range(nb_UD)],axis=0)
    Vis2 = (CF/F_tot)**2
   
    # plt.figure()
    # plt.scatter(q_UD_fitting,Vis2,s=2)
    
    if REG_method == 'TV':

        f_prior = total_variation(HP, I_tot_norm, np.array(res_UD)/2)    

    elif REG_method == 'QS':

        f_prior = quad_smoothness(HP, I_tot_norm)    
        
    residuals = (Vis2-V2_MATISSE)/V2_MATISSE_ERR 
    chi2_tot   = np.sum(residuals**2) + f_prior


    return chi2_tot

def fitting_function(i, wavel_UD, wavel_TOT, HP, diam_UD_tmp, \
                     params2_LM, params2_N, \
                         q_UD_HR_tmp, q_TOT, \
                             intensity_profile_LM, intensity_profile_N, \
                                 flux_ratio_LM, flux_ratio_N, \
                                     V2_TOT, V2_ERR, V_model, \
                                         PATH_OUTPUT_VIS, PATH_OUTPUT_HIST, \
                                             name_var_UD, name_var_F):
    
    if wavel_UD[i]<6:
        diam_UD = diam_UD_tmp[0]
        params2 = params2_LM
        q_UD_HR = q_UD_HR_tmp[0]
        intensity_init = intensity_profile_LM
        flux_init = flux_ratio_LM
        
    else : 
        diam_UD = diam_UD_tmp[1]
        params2 = params2_N
        q_UD_HR = q_UD_HR_tmp[1]
        intensity_init = intensity_profile_N            
        flux_init = flux_ratio_N

    nb_UD = np.shape(diam_UD)[0]        
    
    cond= wavel_TOT == wavel_UD[i]*1E-6
    q_UD_fitting      = np.array(q_TOT[cond])  
    V2_MATISSE     = np.array(V2_TOT[cond])
    V2_MATISSE_ERR = np.array(V2_ERR[cond])     
    
    if nb_UD != 1:
                    

        res2 = minim(short_model, params2, args=(q_UD_fitting, V2_MATISSE, V2_MATISSE_ERR, HP, V_model[i], REG_method, intensity_init, flux_init), method='COBYLA',max_nfev=max_iterations, tol = tolerance) # 5E-4 OK
    
        res_UD = [res2.params['diam_UD'+str(k)].value for k in range(nb_UD)]
        res_UD_bef = [res2.params['diam_UD_bef'+str(k)].value for k in range(nb_UD)]
        res_F_tmp = [res2.params['flux_ftot'+str(k)].value for k in range(nb_UD)]
        
        # This function has been already included in the computation of the chi2 but have to be also execute in the output to take into consideration the different changes made on the flux values.
        # Here we say that wa cannot have an intensity higher than the central uniform disk. Which is in the continuity of our initial condition.
        
        res_F, I_tot_norm = check_intensity(res_F_tmp, res_UD, res_UD_bef, intensity_init, flux_init)


        # Indeed, it is an easy task to compute I_norm (the normalized intensity profil with respect to the central uniform disk value) from the flux since we control the geometry of our model.
        # However in order to gain some computation time and since we have already computed I_norm, we set both as inputs instead of only one of the two.

        # Then based on this new flux ratio for each rings, we compute the coherent flux in order to compute the visibilities

        CF_fitting = np.sum([V_model[i][k]*res_F[k]\
                    for k in range(nb_UD)], axis=0)

                
        F_tot = np.sum([res_F[k] for k in range(nb_UD)],axis=0)
        V2_fitting = (CF_fitting/F_tot)**2
    
        #  The FIT is now over for this fiducial wavlenght, the fitting report is displayed and saved
        
        # resultat=np.append(resultat, res2)
        # print(fit_report(resultat[i]))

        # Here based on the final fitting parameters for this model at this wl, for display purpose, we want to increase the resolution of the model plots
        # Then we compute a last time the visibility curves.

        CF = np.sum([V_ring(q_UD_HR,res_UD_bef[k],res_UD[k])*res_F[k]\
                    for k in range(nb_UD)],axis=0)

        V2_model_HR = (CF/F_tot)**2


        # PLOT of the V2 with higher resolution for the model displayed in linear scale            

        plot_V2_obs_model(q_UD_fitting,V2_MATISSE, q_UD_HR, V2_model_HR, wavel_UD[i], PLOT = False,\
                          SAVE_OUTPUT = PATH_OUTPUT_VIS, logy = True, y_obs_err = V2_MATISSE_ERR,\
                              xlim_min = 1E6, xlim_max = 5E7, ylim_min = 1E-7, ylim_max = 1E0)           

        # PLOT of the V2 with higher resolution for the model displayed in logarithmic scale            
    
        plot_V2_obs_model(q_UD_fitting,V2_MATISSE, q_UD_HR, V2_model_HR, wavel_UD[i], PLOT = False,\
                          SAVE_OUTPUT = PATH_OUTPUT_VIS, y_obs_err = V2_MATISSE_ERR,\
                              xlim_min = 0, xlim_max = 5E7, ylim_min = -0.1, ylim_max = 1.05)           

   
        # Here we set an histogram of the residuals weighted by the error in order to check if we have a Gaussian like distribution to quantify the goodness of fit in addition to the chi2 value
        # However knowing that we use a Bayesian approach the more we increase the prior impact on the fit, the more we will observe an expected deviation from the Gaussian like distribution
   
        plot_histogram(V2_MATISSE,V2_MATISSE_ERR,V2_fitting, wavel_UD[i], PLOT = False, SAVE_OUTPUT = PATH_OUTPUT_HIST)
        
        
        # Now we manually compute the chi2 value and chi2_reduced of the best model fitting at this wavelength taking into account the regularization 
        
        # CHI2:

        N = len(V2_MATISSE)
        chi2_tmp = np.sum(((V2_MATISSE-V2_fitting)/V2_MATISSE_ERR)**2)
        # chi2[i]=chi2_tmp

        # CHI2_REDUCED:

        if nb_UD>1 :
            D_number = N-nb_UD 
            chi2_red=chi2_tmp/D_number
            
        else:
            None
            
        # Regularization:
            
        if REG_method == 'TV':
    
            f_prior_tmp = total_variation(HP, I_tot_norm, np.array(res_UD)/2)    
    
        elif REG_method == 'QS':
    
            f_prior_tmp = quad_smoothness(HP, I_tot_norm)    
               
        # CHI2_TOT
        
        chi2_tot_tmp = chi2_tmp + f_prior_tmp


    # In the case of a simple uniform disk model no fitting at the moment only display

    else:

        params2 = Parameters()        
        params2.add(name_var_UD[0],value = diam_UD[0], vary= False)
        params2.add(name_var_F[0],value = 1, min = 0 , max = 1, vary= False)       
        
        CF = V_uniform(q_UD_HR,params2[name_var_UD[0]].value)*params2[name_var_F[0]].value
        CF_fitting = V_uniform(q_UD_fitting,params2[name_var_UD[0]].value)*params2[name_var_F[0]].value
        V2_model_HR = (np.array(CF))**2
        V2_fitting = (CF_fitting/F_tot)**2


        plot_V2_obs_model(q_UD_fitting,V2_MATISSE, q_UD_HR, V2_model_HR, wavel_UD[i], PLOT = False,\
                          SAVE_OUTPUT = PATH_OUTPUT_VIS, logy = True, y_obs_err = V2_MATISSE_ERR,\
                              xlim_min = 1E6, xlim_max = 5E7, ylim_min = 1E-7, ylim_max = 1E0)           


        plot_V2_obs_model(q_UD_fitting,V2_MATISSE, q_UD_HR, V2_model_HR, wavel_UD[i], PLOT = False,\
                          SAVE_OUTPUT = PATH_OUTPUT_VIS, y_obs_err = V2_MATISSE_ERR,\
                              xlim_min = 0, xlim_max = 5E7, ylim_min = -0.1, ylim_max = 1.05)           



    # INTENSITY PROFILE
    
    res_F_mod = np.array(res_F)/np.sum(res_F) 

    if wavel_UD[i]<6:
        
        res_F = np.round(res_F_mod,6)
        # table_res_fit_LM.add_row([wavel_UD[i]]+ np.round(res_F_mod,6).tolist())
        I_tot = I_tot_norm


    else :
        res_F = np.round(res_F_mod,6)
        I_tot = I_tot_norm

    print(fit_report(res2))

    print('WORK ON %.3f OVER'%wavel_UD[i])

    
    return  wavel_UD[i], res2, chi2_tmp, chi2_red, f_prior_tmp, chi2_tot_tmp, res_F, I_tot

def UD_modeling(wavel_UD,wavel_TOT, q_TOT, V2_TOT, V2_ERR,\
                         q_UD_HR_tmp, diam_UD_tmp, diam_UD_bef_tmp, flux_ratio_LM, flux_ratio_N, intensity_profile_LM, intensity_profile_N, flux_LM_max_tmp, flux_N_max_tmp,\
                                     HP, V_model,\
                                         PATH_OUTPUT, PATH_OUTPUT_VIS, PATH_OUTPUT_INT, PATH_OUTPUT_HIST, PATH_OUTPUT_FIT_RES, PLOT):
         

    I_tot_LM = np.zeros((sum(i < 6 for i in wavel_UD),len(diam_UD_tmp[0])))
    I_tot_N  = np.zeros((sum(i > 6 for i in wavel_UD),len(diam_UD_tmp[1])))
    
    resultat=[]

    # INITIALIZATION OF THE RESULTS TABLE OF THE FITTING FLUX FOR EACH COMPONENTS OF THE MODEL: DISPLAY ==> FLUX vs RADIUS 
    
    table_res_fit_LM = PrettyTable()
    table_res_fit_LM.field_names = ["WAVEL"] + [str((np.array(diam_UD_tmp)/2)[0][i]) for i in range(len(diam_UD_tmp[0]))]

    table_res_fit_N = PrettyTable()
    table_res_fit_N.field_names = ["WAVEL"] + [str((np.array(diam_UD_tmp)/2)[1][i]) for i in range(len(diam_UD_tmp[1]))]
    

    # INITIALIZATION OF THE MINIMIZATION ROUTINE PARAMETER'S 

    # FOR LM-BAND

    name_var_UD = ['diam_UD'+str(i) for i in range(np.shape(diam_UD_tmp[0])[0])]
    name_var_UD_bef = ['diam_UD_bef'+str(i) for i in range(np.shape(diam_UD_tmp[0])[0])]
    name_var_F = ['flux_ftot'+str(i) for i in range(np.shape(diam_UD_tmp[0])[0])]

    variation_LM = np.ones(len(diam_UD_tmp[0])).astype('bool')
    variation_LM[0] = False

    params2_LM = Parameters()        
    [params2_LM.add(name_var_UD[k],value = diam_UD_tmp[0][k], vary= False) for k in range(np.shape(diam_UD_tmp[0])[0])]
    [params2_LM.add(name_var_UD_bef[k],value = diam_UD_bef_tmp[0][k], vary= False) for k in range(np.shape(diam_UD_tmp[0])[0])]
    [params2_LM.add(name_var_F[k],value = flux_ratio_LM[k], min = 0 , max = flux_LM_max_tmp[k], vary= variation_LM[k]) for k in range(np.shape(diam_UD_tmp[0])[0])]      

    # Indeed the more we get further from the star, the more the rings have low probability to have higher intensity than the central uniform disk using this maximum boundary
    # This will also hardly depend on the initial intensity profile that we use. At the moment only 1/r^alpha are available thus this is impossible 
    # to have high intensity rings far from the star.


    # FOR N-BAND

    name_var_UD = ['diam_UD'+str(i) for i in range(np.shape(diam_UD_tmp[1])[0])]
    name_var_UD_bef = ['diam_UD_bef'+str(i) for i in range(np.shape(diam_UD_tmp[1])[0])]
    name_var_F = ['flux_ftot'+str(i) for i in range(np.shape(diam_UD_tmp[1])[0])]

    variation_N = np.ones(len(diam_UD_tmp[1])).astype('bool')
    variation_N[0] = False

    params2_N = Parameters()        
    [params2_N.add(name_var_UD[k],value = diam_UD_tmp[1][k], vary= False) for k in range(np.shape(diam_UD_tmp[1])[0])]
    [params2_N.add(name_var_UD_bef[k],value = diam_UD_bef_tmp[1][k], vary= False) for k in range(np.shape(diam_UD_tmp[1])[0])]
    [params2_N.add(name_var_F[k],value = flux_ratio_N[k], min = 0 , max = flux_N_max_tmp[k], vary= variation_N[k]) for k in range(np.shape(diam_UD_tmp[1])[0])]      

    
    # OK LET'S GO FIT: DON'T STOP ME NOW !

    
    pool = mp.Pool(processes=nprocs)    
        
    result_parallel = pool.starmap(fitting_function, [(i, wavel_UD, wavel_TOT, HP, diam_UD_tmp, \
                         params2_LM, params2_N, \
                             q_UD_HR_tmp, q_TOT, \
                                 intensity_profile_LM, intensity_profile_N, \
                                     flux_ratio_LM, flux_ratio_N, \
                                         V2_TOT, V2_ERR, V_model, \
                                             PATH_OUTPUT_VIS, PATH_OUTPUT_HIST, \
                                                 name_var_UD, name_var_F) for i in range(len(wavel_UD))])
    
    pool.close()
    pool.join()



    wavel_UD, resultat, chi2, chi2_red, f_prior, chi2_tot, res_F, I_tot = np.array(result_parallel).T

    # TABLE OF RESULTS

    i_LM = 0
    i_N  = 0

    for ind in range(len(wavel_UD)):        
        
        if wavel_UD[ind]<6:
            table_res_fit_LM.add_row([wavel_UD[ind]] + np.round(res_F[ind],6).tolist())
            I_tot_LM[i_LM] = I_tot[ind]
            i_LM = i_LM+1
    
    
        else :
            table_res_fit_N.add_row([wavel_UD[ind]]+ np.round(res_F[ind],6).tolist())
            I_tot_N[i_N] = I_tot[ind]
            i_N = i_N+1


    # TABLE OF FLUX

    with open(PATH_OUTPUT_FIT_RES+'fit_flux_ratio_LM.dat', 'w') as f: f.write(str(table_res_fit_LM))
    with open(PATH_OUTPUT_FIT_RES+'fit_flux_ratio_N.dat', 'w') as f: f.write(str(table_res_fit_N))

    # TABLE OF INTENSITY
    wavel_UD = np.array(wavel_UD)
    data_intensity = PrettyTable()
    column_names_I = ["rho_ext"]+[str(np.round(wavel_UD[wavel_UD<6][i],5)) for i in range(len(I_tot_LM))]
    data_intensity.add_column(column_names_I[0], (np.array(diam_UD_tmp[0])/2).tolist())
    for i in range(len(I_tot_LM)):
        data_intensity.add_column(column_names_I[i+1],I_tot_LM[i])
    with open(PATH_OUTPUT_FIT_RES+'intensity_LM.dat', 'w') as f: f.write(str(data_intensity))

    data_intensity = PrettyTable()
    column_names_I = ["rho_ext"]+[str(np.round(wavel_UD[wavel_UD>=6][i],5)) for i in range(len(I_tot_N))]
    data_intensity.add_column(column_names_I[0], (np.array(diam_UD_tmp[1])/2).tolist())
    for i in range(len(I_tot_N)):
        data_intensity.add_column(column_names_I[i+1],I_tot_N[i])
    with open(PATH_OUTPUT_FIT_RES+'intensity_N.dat', 'w') as f: f.write(str(data_intensity))


    # ALL LMFIT REPORT

    with open(PATH_OUTPUT_FIT_RES+ 'resultat_fitting_rings'+'.txt', "w+") as att_file:
        for ind in range(len(wavel_UD)):
            att_file.write(str(wavel_UD[ind])+ ' µm' + '\n')
            att_file.write(fit_report(resultat[ind]) + '\n')

    # ALL CHI2 REPORT
        
    data_chi2 = PrettyTable()
    data_chi2.field_names = ['Wavel', 'chi2', 'chi2_red', 'f_prior_'+REG_method, 'chi2_tot']
    for k in range(len(chi2)):
        data_chi2.add_row([wavel_UD[k]]+ [chi2[k]] + [chi2_red[k]] + [f_prior[k]] + [chi2_tot[k]])
    with open(PATH_OUTPUT_FIT_RES+'chi2.dat', 'w') as f: f.write(str(data_chi2))
    
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









