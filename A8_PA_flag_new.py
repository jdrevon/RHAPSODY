#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 14:55:09 2020

@author: jdrevon
"""
import numpy as np

def flag_PA_new(x):  

    flag = np.ones(len(x))*np.nan
    bornes_test = np.linspace(-180,180,18+1)
    index_0_test = int(np.where(bornes_test==0)[0][0])
    for k in range(len(x)):        
        for i in range(index_0_test):
            cond1 = bornes_test[i] <= x[k] < bornes_test[i+1]
            cond2 = bornes_test[i]+180 <= x[k] < bornes_test[i+1]+180
            if np.logical_or(cond1,cond2)==True :
                flag[k]=int(i) 
    return flag
