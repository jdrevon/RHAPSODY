# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:16:19 2022

@author: jdrevon
"""

import numpy as np

def uniform_disk(diam_list,diam):

    I_ratio = np.zeros(len(diam_list))
    I_ratio[diam_list/2<diam/2]= 1
    
    return I_ratio