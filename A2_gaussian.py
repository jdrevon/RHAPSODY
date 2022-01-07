#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 18:53:30 2020

@author: Julien
"""

import numpy as np

def gaussian(central_value,sigma,A,x):
            
    gaussian_value = A*np.exp(-1/2*((x-central_value)/sigma)**2)
    
    return gaussian_value

def gaussian_exp(x,central_value,sigma,A):
            
    gaussian_value = A*np.exp(-1/2*((x-central_value)/sigma)**2)
    
    return np.exp(-gaussian_value)


# /(sigma*np.sqrt(2*np.pi))