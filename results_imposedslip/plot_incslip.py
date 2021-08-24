#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:25:04 2020
Script to read output from BEM folding models and plot incremental slip
@author: rishav
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.collections import PolyCollection
plotrange = 10
#%%
numstart = 89#term
numend = numstart
for num in range(numstart,numend+1):
    print(num)
    ofol = "modelrun_" + str(num)

    nfil=0
    flist = []
    totalslip = 1e-6

    for file in os.listdir(ofol):
        if file.startswith("data"):
            print(os.path.join(ofol, file))
            flist.insert(nfil,ofol+file)
            nfil=nfil+1
            
            metafile = ofol + "/" + "inputs.dat"
            datafile = ofol + "/" + file
            
            
            fileout = ofol + "/incslip_" + file[0:len(file)-4] + ".jpg"
            print(fileout)
            
            metadata=pd.read_table(metafile,header=None)
            datatable=pd.read_table(datafile,header=None)

            
            
            data = datatable.values

            dip = np.arctan2(-(data[:,3]-data[:,1]),data[:,2]-data[:,0])
            pdip = dip+np.pi/2
            w = 0.02
            
            vertices = []
            for i in range(0,len(data[:,0])-1):
                vertices.append([(data[i,0]+w*np.cos(pdip[i]),data[i,1]-w*np.sin(pdip[i])),
                                 (data[i,2]+w*np.cos(pdip[i]),data[i,3]-w*np.sin(pdip[i])),
                                 (data[i,2]-w*np.cos(pdip[i]),data[i,3]+w*np.sin(pdip[i])),
                                 (data[i,0]-w*np.cos(pdip[i]),data[i,1]+w*np.sin(pdip[i]))])
            
            # Plot Slip on fault (ratio wrt total slip)
            z = 100*data[:,8]/np.amax(data[:,8])

            fig = plt.figure(1,figsize=(30,4))
            fig.clf()
            ax = fig.add_subplot(1,1,1)
            im = PolyCollection(vertices,array=z,cmap='Spectral_r')
            ax.add_collection(im)
            plt.axis('tight')
            ax.set_aspect('equal','box')
            cb=plt.colorbar(im,ax=ax)
            cb.ax.set_ylabel('RL ----- Slip% ----- LL')
            ax.set_xlabel('X(km)')
            ax.set_ylabel('Z(km)')
            im.set_clim(-plotrange,plotrange)
            # im.set_clim(0,plotrange)
            plt.xlim(-12,12)
            plt.ylim(-3,0)
            #plt.title(titlelabel)
            #plt.show()

            plt.savefig(fileout)

            
            
