# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:47:46 2024

@author: at924
"""

import numpy as np
import matplotlib.pyplot as plt


fpath='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/Random/'

fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Random_DEL/'


cases=['_5ms','_10ms','_15ms','_20ms','_25ms']

cases2=['_5ms','_10ms','_15ms','_20ms','_25ms']

cut=12000
time=np.linspace(0, 720-cut*0.01,72002-cut)
Amax=np.zeros(shape=(7,10))
Amax2=np.zeros(shape=(7,10))

for i in range(0,len(cases)):
    for j in range(1,11):
        data=np.genfromtxt(fpath + cases[i] +  '/aTOP' +cases[i] +'_SET' + str(j) + '.txt')
        #plt.figure(i)
        #plt.plot(time,data[cut:,1],label=str(j))
        Amax[i,j-1]=np.max(np.abs(data[cut:,1]))
        
for i in range(0,len(cases2)):
    for j in range(1,11):
        data=np.genfromtxt(fpath2 + cases2[i] +  '/aTOP' + '_SET' + str(j) +cases2[i] + '_DEL.txt')
        #plt.figure(200+i)
        #plt.plot(time,data[cut:-1],label=str(j))
        Amax2[i,j-1]=np.max(np.abs(data[cut:]))
        
        

plt.figure(100)
xs=[5.,10.,15.,20.,25.]
[plt.plot([xs[i]]*10,Amax[i,:],'*',color='black') for i in range(0,len(xs))]

[plt.plot([xs[i]]*10,Amax2[i,:],'+',color='red') for i in range(0,len(xs))]








        