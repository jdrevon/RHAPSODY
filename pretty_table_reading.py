#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 13:24:33 2021

@author: jdrevon
"""
import numpy as np

def READ_PRETTY_TABLE(PATH,HEADER_LINES):
    
    # PATH ==> Location of the 
    
    result=[]
    with open(PATH) as f:
        s = [line.rstrip('\n') for line in f]
        for line in s:
            splitdata = line.split("|")
            if len(splitdata) == 1:
                continue  # skip lines with no separators
            linedata = []
            for field in splitdata:
                field = field.strip()
                if field:
                    linedata.append(field)
            try:
                result.append(list(map(float, linedata)))
            except:
                result.append(linedata)

#    HEADER_LINES = 1
    
    header = np.array(result[:HEADER_LINES][0])
    data = np.array(result[HEADER_LINES:])
    
    return header,data

