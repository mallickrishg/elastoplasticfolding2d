#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:25:04 2020
Script to read output from BEM folding models and plot it
@author: rishav
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.collections import PolyCollection
plotrange = 10
#%%
numstart = 74
numend = 85#numstart
for num in range(numstart,numend+1):
    print(num)
    ofol = "modelrun_" + str(num)

    metafile = ofol + "/" + "inputs.dat"
    metadata=pd.read_table(metafile,header=None)    
            

    #%% make figure for input parameters
    collabel = ("G (GPa)","nu","mu_l","mu_f","Cohesion (MPa)","Sxx_i","Sxz_i","Szz_i","lambda","Totslip (km)","nslips")
            
    plt.figure(3,figsize=(35,5))
    plt.clf()
    tblvals = metadata.values
    tblvals[0,0] = tblvals[0,0]/1e9
    tblvals[0,4] = tblvals[0,4]/1e6
    tblvals[0,5] = tblvals[0,5]/1e6

    tbl = plt.table(cellText = tblvals,colLabels=collabel,loc='center')
    tbl.set_fontsize(20)
    tbl.scale(1,3)
    plt.text(0.5,0.7,'Parameters',size = 20)
    #plt.show()
    plt.savefig(ofol + "/inputdata.jpg")

