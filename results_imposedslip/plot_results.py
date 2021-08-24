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
numstart = 20
numend = 27#numstart
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
            
            
            fileout = ofol + "/figure_" + file[0:len(file)-4] + ".jpg"
            print(fileout)
            
            metadata=pd.read_table(metafile,header=None)
            datatable=pd.read_table(datafile,header=None)

            
            
            data = datatable.values
            #totalslip = metadata.values[0,9]
            totalslip = data[1,5] 
            #print(totalslip)

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
            z = 100*data[:,4]/(1e3*totalslip)
            #z = data[:,4]
            #plotrange = 100#0.1*np.max(z)

            fig = plt.figure(1,figsize=(30,4))
            fig.clf()
            ax = fig.add_subplot(1,1,1)
            im = PolyCollection(vertices,array=z,cmap='jet')
            ax.add_collection(im)
            plt.axis('tight')
            ax.set_aspect('equal','box')
            cb=plt.colorbar(im,ax=ax)
            cb.ax.set_ylabel('RL ----- Slip% ----- LL')
            ax.set_xlabel('X(km)')
            ax.set_ylabel('Z(km)')
            im.set_clim(-plotrange,plotrange)
            plt.xlim(-12,12)
            #plt.ylim(-3,0)
            #plt.title(titlelabel)
            #plt.show()

            plt.savefig(fileout)

            
            # Plot Stress on fault (ratio wrt yield)
            z = data[:,6]/data[:,7]
            fig = plt.figure(2,figsize=(30,4))
            fig.clf()
            ax = fig.add_subplot(1,1,1)
            im = PolyCollection(vertices,array=z,cmap='jet')
            ax.add_collection(im)
            plt.axis('tight')
            ax.set_aspect('equal','box')
            cb=plt.colorbar(im,ax=ax)
            cb.ax.set_ylabel('tau/tau_y')
            ax.set_xlabel('X(km)')
            ax.set_ylabel('Z(km)')
            im.set_clim(-1.5,1.5)
            plt.xlim(-12,12)
            taufileout = ofol + "/taufigure_" + file[0:len(file)-4] + ".jpg"
            plt.savefig(taufileout)


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







