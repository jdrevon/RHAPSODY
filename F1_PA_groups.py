# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:00:00 2021

@author: jdrevon
"""

import numpy as np
from RHAPSODY_init import NBR_GRPS_PA

def PA_groups(x):
    
    '''    
    Input: I
    
    x : array of number between -180° and 180°
    
    Output: O
    
    PA_groups : array contening the number associated to the PA groups (see terminal for more information)
    
    '''
    
    PA_groups  = np.full(len(x),999)
    PA_boundaries = np.linspace(0,180,NBR_GRPS_PA+1) #Decompose the interval in NBR_GRPS_PA independant groups
    
    test = ["Group %i : %i-%i" %(i,PA_boundaries[i],PA_boundaries[i+1]) for i in range(NBR_GRPS_PA)]
    print("PA Groups:", *test, sep="\n")
    
    cond_0 = (x==0)
    PA_groups[cond_0] = 0
    
    for k in range(NBR_GRPS_PA):
        
        cond = np.logical_and(np.abs(x)>PA_boundaries[k],np.abs(x)<=PA_boundaries[k+1])
        PA_groups[cond] = int(k)
        
    PA_groups = PA_groups.astype('int')

    return PA_groups    