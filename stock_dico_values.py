#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 09:07:41 2021

@author: jdrevon
"""

import numpy.ma as ma

def stock_V2_from_dico(OIFITS_TOT_LM, OIFITS_TOT_N):
    
    V2_TOT = []
    V2_ERR = []
    q_TOT  = []
    wavel_TOT = []
    U_TOT  = []
    V_TOT  = []
     
    for i in range(len(OIFITS_TOT_LM)):
        
        flag           = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['FLAG'])
        q_tmp          = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['BASELINE'].astype('float')/OIFITS_TOT_LM[i]['VIS2']['WAVEL'].astype('float'))
        wavel_tmp      = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['WAVEL'].astype('float'))
        U              = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['U'].astype('float'))
        V              = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['V'].astype('float'))
        V2_tmp         = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['VIS2'].astype('float'))
        V2_err_tmp     = ma.concatenate(OIFITS_TOT_LM[i]['VIS2']['VIS2_ERR'].astype('float'))

        V2_TOT.append(ma.masked_array(V2_tmp,mask=flag))
        V2_ERR.append(ma.masked_array(V2_err_tmp,mask=flag))
        q_TOT.append(ma.masked_array(q_tmp,mask=flag))
        wavel_TOT.append(ma.masked_array(wavel_tmp,mask=flag))    
        U_TOT.append(ma.masked_array(U,mask=flag))        
        V_TOT.append(ma.masked_array(V,mask=flag))     
        
    for k in range(len(OIFITS_TOT_N)):
        
        flag           = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['FLAG'])
        q_tmp          = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['BASELINE'].astype('float')/OIFITS_TOT_N[k]['VIS2']['WAVEL'].astype('float'))
        wavel_tmp      = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['WAVEL'].astype('float'))
        U              = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['U'].astype('float'))
        V              = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['V'].astype('float'))
        V2_tmp         = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['VIS2'].astype('float'))
        V2_err_tmp     = ma.concatenate(OIFITS_TOT_N[k]['VIS2']['VIS2_ERR'].astype('float'))


        
        
        V2_TOT.append(ma.masked_array(V2_tmp,mask=flag))
        V2_ERR.append(ma.masked_array(V2_err_tmp,mask=flag))
        q_TOT.append(ma.masked_array(q_tmp,mask=flag))
        wavel_TOT.append(ma.masked_array(wavel_tmp,mask=flag))           
        
        
    V2_TOT = ma.concatenate(V2_TOT)
    V2_ERR = ma.concatenate(V2_ERR)
    q_TOT  = ma.concatenate(q_TOT)
    wavel_TOT = ma.concatenate(wavel_TOT) 
    
    return wavel_TOT.compressed(), q_TOT.compressed(), V2_TOT.compressed(), V2_ERR.compressed()


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