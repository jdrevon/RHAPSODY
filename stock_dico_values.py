#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 09:07:41 2021

@author: jdrevon
"""

import numpy.ma as ma
from rhapsody_init import DATA_band_name
import numpy as np

def stock_V2_from_dico(OIFITS_TOT):
    
    V2_TOT = np.zeros(len(OIFITS_TOT),dtype=object)
    wavel_TOT = np.zeros(len(OIFITS_TOT),dtype=object)
    q_TOT = np.zeros(len(OIFITS_TOT),dtype=object)
    V2_ERR = np.zeros(len(OIFITS_TOT),dtype=object)
    V_TOT = np.zeros(len(OIFITS_TOT),dtype=object)
    U_TOT = np.zeros(len(OIFITS_TOT),dtype=object)

    
    for k in range(len(OIFITS_TOT)):

        V2_TOT_tmp = []
        V2_ERR_tmp = []
        q_TOT_tmp  = []
        wavel_TOT_tmp = []
        U_TOT_tmp  = []
        V_TOT_tmp  = []
    
    
        for i in range(len(OIFITS_TOT[k])):
            
            try:
                flag           = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['FLAG'])
                q_tmp          = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['BASELINE'].astype('float')/OIFITS_TOT[k][i]['VIS2']['WAVEL'].astype('float'))
                wavel_tmp      = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['WAVEL'].astype('float'))
                U              = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['U'].astype('float'))
                V              = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['V'].astype('float'))
                V2_tmp         = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['VIS2'].astype('float'))
                V2_err_tmp     = ma.concatenate(OIFITS_TOT[k][i]['VIS2']['VIS2_ERR'].astype('float'))

            except: # In the cas only 1 wavelenght is present
                
                print('INFO: ONLY 1 wavelenght is detected in the file')
                
                flag           = OIFITS_TOT[k][i]['VIS2']['FLAG']
                q_tmp          = OIFITS_TOT[k][i]['VIS2']['BASELINE'].astype('float')/OIFITS_TOT[k][i]['VIS2']['WAVEL'].astype('float')
                wavel_tmp      = OIFITS_TOT[k][i]['VIS2']['WAVEL'].astype('float')
                U              = OIFITS_TOT[k][i]['VIS2']['U'].astype('float')
                V              = OIFITS_TOT[k][i]['VIS2']['V'].astype('float')
                V2_tmp         = OIFITS_TOT[k][i]['VIS2']['VIS2'].astype('float')
                V2_err_tmp     = OIFITS_TOT[k][i]['VIS2']['VIS2_ERR'].astype('float')
    
            V2_TOT_tmp.append(ma.masked_array(V2_tmp,mask=flag))
            V2_ERR_tmp.append(ma.masked_array(V2_err_tmp,mask=flag))
            q_TOT_tmp.append(ma.masked_array(q_tmp,mask=flag))
            wavel_TOT_tmp.append(ma.masked_array(wavel_tmp,mask=flag))    
            U_TOT_tmp.append(ma.masked_array(U,mask=flag))        
            V_TOT_tmp.append(ma.masked_array(V,mask=flag))     
        
        U_TOT_tmp = ma.concatenate(U_TOT_tmp)        
        V_TOT_tmp = ma.concatenate(V_TOT_tmp)     

        V2_TOT_tmp = ma.concatenate(V2_TOT_tmp)
        V2_ERR_tmp = ma.concatenate(V2_ERR_tmp)
        q_TOT_tmp  = ma.concatenate(q_TOT_tmp)
        wavel_TOT_tmp = ma.concatenate(wavel_TOT_tmp) 


        V2_TOT[k] = V2_TOT_tmp.compressed()
        q_TOT[k] = q_TOT_tmp.compressed()
        V2_ERR[k] = V2_ERR_tmp.compressed()
        wavel_TOT[k] = wavel_TOT_tmp.compressed()
        U_TOT [k]  = U_TOT_tmp.compressed()
        V_TOT [k]  = V_TOT_tmp.compressed()
    
    return wavel_TOT, q_TOT, V2_TOT, V2_ERR, U_TOT, V_TOT
        


def stock_TP_from_dico(OIFITS_TOT_LM, OIFITS_TOT_N):

    wavel_T3_TOT, B1_TOT,B2_TOT,B3_TOT,T3_TOT,T3_err_TOT, U1_T_TOT,U2_T_TOT,U3_T_TOT,V1_T_TOT,V2_T_TOT,V3_T_TOT =[],[],[],[],[],[],[],[],[],[],[],[]
     
    for i in range(len(OIFITS_TOT_LM)):
        
        flag_T3        = ma.concatenate(OIFITS_TOT_LM[i]['T3']['FLAG'])
        wavel_T3       = ma.concatenate(OIFITS_TOT_LM[i]['T3']['WAVEL'].astype('float'))
   
        U1             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['U1'].astype('float'))
        U2             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['U2'].astype('float'))
        U3             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['U3'].astype('float'))

        V1             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['V1'].astype('float'))
        V2             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['V2'].astype('float'))
        V3             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['V3'].astype('float'))
        
   
        B1             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['B1'].astype('float'))
        B2             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['B2'].astype('float'))
        B3             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['B3'].astype('float'))
        T3             = ma.concatenate(OIFITS_TOT_LM[i]['T3']['T3'].astype('float'))
        T3_err         = ma.concatenate(OIFITS_TOT_LM[i]['T3']['T3_ERR'].astype('float'))


        # print(np.shape(U),np.shape(q_tmp))  
        
        B1_TOT.append(ma.masked_array(B1,mask=flag_T3))   
        B2_TOT.append(ma.masked_array(B2,mask=flag_T3))   
        B3_TOT.append(ma.masked_array(B3,mask=flag_T3))   

        U1_T_TOT.append(ma.masked_array(U1,mask=flag_T3))   
        U2_T_TOT.append(ma.masked_array(U2,mask=flag_T3))   
        U3_T_TOT.append(ma.masked_array(U3,mask=flag_T3))   

        V1_T_TOT.append(ma.masked_array(V1,mask=flag_T3))   
        V2_T_TOT.append(ma.masked_array(V2,mask=flag_T3))   
        V3_T_TOT.append(ma.masked_array(V3,mask=flag_T3))   


        T3_TOT.append(ma.masked_array(T3,mask=flag_T3))   
        T3_err_TOT.append(ma.masked_array(T3_err,mask=flag_T3))   
        wavel_T3_TOT.append(ma.masked_array(wavel_T3,mask=flag_T3))    

        
    for k in range(len(OIFITS_TOT_N)):

        flag_T3        = ma.concatenate(OIFITS_TOT_N[k]['T3']['FLAG'])
        wavel_T3      = ma.concatenate(OIFITS_TOT_N[k]['T3']['WAVEL'].astype('float'))
        
        B1             = ma.concatenate(OIFITS_TOT_N[k]['T3']['B1'].astype('float'))
        B2             = ma.concatenate(OIFITS_TOT_N[k]['T3']['B2'].astype('float'))
        B3             = ma.concatenate(OIFITS_TOT_N[k]['T3']['B3'].astype('float'))
        
        U1             = ma.concatenate(OIFITS_TOT_N[k]['T3']['U1'].astype('float'))
        U2             = ma.concatenate(OIFITS_TOT_N[k]['T3']['U2'].astype('float'))
        U3             = ma.concatenate(OIFITS_TOT_N[k]['T3']['U3'].astype('float'))

        V1             = ma.concatenate(OIFITS_TOT_N[k]['T3']['V1'].astype('float'))
        V2             = ma.concatenate(OIFITS_TOT_N[k]['T3']['V2'].astype('float'))
        V3             = ma.concatenate(OIFITS_TOT_N[k]['T3']['V3'].astype('float'))

        
        T3             = ma.concatenate(OIFITS_TOT_N[k]['T3']['T3'].astype('float'))
        T3_err         = ma.concatenate(OIFITS_TOT_N[k]['T3']['T3_ERR'].astype('float'))
   
    
        B1_TOT.append(ma.masked_array(B1,mask=flag_T3))   
        B2_TOT.append(ma.masked_array(B2,mask=flag_T3))   
        B3_TOT.append(ma.masked_array(B3,mask=flag_T3)) 

        U1_T_TOT.append(ma.masked_array(U1,mask=flag_T3))   
        U2_T_TOT.append(ma.masked_array(U2,mask=flag_T3))   
        U3_T_TOT.append(ma.masked_array(U3,mask=flag_T3))   

        V1_T_TOT.append(ma.masked_array(V1,mask=flag_T3))   
        V2_T_TOT.append(ma.masked_array(V2,mask=flag_T3))   
        V3_T_TOT.append(ma.masked_array(V3,mask=flag_T3))   

        T3_TOT.append(ma.masked_array(T3,mask=flag_T3))   
        T3_err_TOT.append(ma.masked_array(T3_err,mask=flag_T3))   
        wavel_T3_TOT.append(ma.masked_array(wavel_T3,mask=flag_T3))
        
    B1_TOT = ma.concatenate(B1_TOT)
    B2_TOT = ma.concatenate(B2_TOT)
    B3_TOT =ma.concatenate(B3_TOT)
    T3_TOT =ma.concatenate(T3_TOT)
    T3_err_TOT =ma.concatenate(T3_err_TOT)
    wavel_T3_TOT =ma.concatenate(wavel_T3_TOT)
    
    U1_T_TOT = ma.concatenate(U1_T_TOT)   
    U2_T_TOT = ma.concatenate(U2_T_TOT)   
    U3_T_TOT = ma.concatenate(U3_T_TOT)   

    V1_T_TOT = ma.concatenate(V1_T_TOT)   
    V2_T_TOT = ma.concatenate(V2_T_TOT)   
    V3_T_TOT = ma.concatenate(V3_T_TOT)   

    
    return B1_TOT,B2_TOT,B3_TOT,T3_TOT,T3_err_TOT, wavel_T3_TOT, U1_T_TOT,U2_T_TOT,U3_T_TOT,V1_T_TOT,V2_T_TOT,V3_T_TOT