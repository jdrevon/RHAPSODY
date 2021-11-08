# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 19:19:11 2021

@author: jdrevon
"""
import numpy as np
from brightness_distribution_models import normalized_brightness_UD, normalized_brightness_R
from copy import deepcopy

def check_intensity_old(res_UD,res_UD_bef,res_F_mod,nb_UD):

        rho   = np.array(res_UD)/2

        # I_UD_norm = normalized_brightness_UD(np.array([0]),res_UD[0])/normalized_brightness_UD(np.array([0]),res_UD[0])[0]
        F_UD_0 = res_F_mod[0]
        F_new = deepcopy(res_F_mod)
        I_tot_norm = np.sum([res_F_mod[k]/F_UD_0\
                                         *normalized_brightness_R(rho, res_UD_bef[k], res_UD[k])/normalized_brightness_UD(rho,res_UD[0])[0]\
                                         for k in range(nb_UD)],axis=0)

        # print(I_tot_norm)
        r_prob_index = np.where(I_tot_norm>1)[0]
          
        
        if np.shape(r_prob_index)[0]!= 0:
            # print('rho \n')
            # print(rho)
            # print('I_tot_norm \n')
            # print(I_tot_norm)
            # print('r_prob_index \n')
            # print(r_prob_index)
            # print('res_F_mod \n')
            # print(res_F_mod)
            for i in range(0,len(r_prob_index)):
                index = r_prob_index[i]
                ratio_I = normalized_brightness_UD(np.array([0]),res_UD[0])[0]/normalized_brightness_R(np.array([rho[index]]), res_UD_bef[index], res_UD[index])[0]
                F_new[index] = ratio_I*res_F_mod[0]
                I_tot_norm[index] = np.sum(F_new[index]/F_UD_0/ratio_I,axis=0)
            # print(F_new)
            # print('I_tot_norm_NEW \n')
            # print(I_tot_norm)
            
        return F_new, I_tot_norm


def check_intensity(res_F_mod, res_UD, res_UD_bef, intensity_init, flux_init):

        I_tot_tmp = (res_F_mod*np.array([normalized_brightness_R(res_UD[i]/2, res_UD_bef[i], res_UD[i])[0] for i in range(len(res_UD))]))
        I_tot_norm = I_tot_tmp#/I_tot_tmp[0]

        r_prob_index = np.where(I_tot_norm>1)[0]
                  
        F_new = res_F_mod
                
        if np.shape(r_prob_index)[0]!= 0:
            for i in range(0,len(r_prob_index)):
                index = r_prob_index[i]
                F_new[index] = flux_init[index]
                I_tot_norm[index] =  intensity_init[index]

        #F_new[index] = 1/normalized_brightness_R(res_UD[index]/2, res_UD_bef[index], res_UD[index])[0]
                

        # F_UD_0 = res_F_mod[0]
        # F_new = res_F_mod
        # I_tot_norm = res_F_mod/F_UD_0*intensity_init/intensity_init[0]

        # # print(I_tot_norm)
        # r_prob_index = np.where(I_tot_norm>1)[0]
          
        
        # if np.shape(r_prob_index)[0]!= 0:
        #     # print('rho \n')
        #     # print(rho)
        #     # print('I_tot_norm \n')
        #     # print(I_tot_norm)
        #     # print('r_prob_index \n')
        #     # print(r_prob_index)
        #     # print('res_F_mod \n')
        #     # print(res_F_mod)
        #     for i in range(0,len(r_prob_index)):
        #         index = r_prob_index[i]
        #         ratio_I = intensity_init[0]/intensity_init[index]
        #         F_new[index] = ratio_I*res_F_mod[0]
        #         I_tot_norm[index] = np.sum(F_new[index]/F_UD_0/ratio_I,axis=0)
                
            
        return F_new, I_tot_norm