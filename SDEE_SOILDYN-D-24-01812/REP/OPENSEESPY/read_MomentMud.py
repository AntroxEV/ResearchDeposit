# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:47:46 2024

@author: at924
"""

import numpy as np
import matplotlib.pyplot as plt


fpath='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/Random/'

fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Random/'


cases=['_5ms','_10ms','_12_5ms','_15ms','_17_5ms','_20ms','_25ms']

cases2=['_5ms','_10ms','_15ms','_17_5ms','_20ms','_25ms']

cut=12000
time=np.linspace(0, 720-cut*0.01,72002-cut)
MYmax=np.zeros(shape=(7,10))
MYmax2=np.zeros(shape=(7,10))

for i in range(0,len(cases)):
    for j in range(1,11):
        data=np.genfromtxt(fpath + cases[i] +  '/MYmud' +cases[i] +'_SET' + str(j) + '.txt')
        #plt.figure(i)
        #plt.plot(time,data[cut:,1],label=str(j))
        MYmax[i,j-1]=np.max(np.abs(data[cut:,1]))
        
for i in range(0,len(cases2)):
    for j in range(1,11):
        data=np.genfromtxt(fpath2 + cases2[i] +  '/MYmud' + '_SET' + str(j) +cases2[i] + '_OPT.txt')
        #plt.figure(200+i)
        #plt.plot(time,data[cut:-1],label=str(j))
        MYmax2[i,j-1]=np.max(np.abs(data[cut:]))
        
        

plt.figure(100)
xs=[5.,10.,12.5,15.,17.5,20.,25.]
plt.plot([xs[0]]*10,MYmax[0,:],'*',color='black')
plt.plot([xs[1]]*10,MYmax[1,:],'*',color='black')
plt.plot([xs[2]]*10,MYmax[2,:],'*',color='black')
plt.plot([xs[3]]*10,MYmax[3,:],'*',color='black')
plt.plot([xs[4]]*10,MYmax[4,:],'*',color='black')
plt.plot([xs[5]]*10,MYmax[5,:],'*',color='black')
plt.plot([xs[6]]*10,MYmax[6,:],'*',color='black')

xs=[5.,10.,15.,17.5,20.,25.]
plt.plot([xs[0]]*10,MYmax2[0,:],'+',color='red')
plt.plot([xs[1]]*10,MYmax2[1,:],'+',color='red')
plt.plot([xs[2]]*10,MYmax2[2,:],'+',color='red')
plt.plot([xs[3]]*10,MYmax2[3,:],'+',color='red')
plt.plot([xs[4]]*10,MYmax2[4,:],'+',color='red')




for i in range(0,1):
    for j in range(1,11):
        plt.figure()
        data=np.genfromtxt(fpath + cases2[i] +  '/MYmud' +cases2[i] +'_SET' + str(j) + '.txt')
        data2=np.genfromtxt(fpath2 + cases2[i] +  '/MYmud' + '_SET' + str(j) +cases2[i] + '_OPT.txt')
        plt.plot(time,data[cut:,1],'k',label='NL case' + cases2[i] +'SET {:.0f}'.format(j))
        plt.plot(time,data2[cut:-1],'--r',label='ADamp err={:.2%}'.format(1-np.divide(MYmax2[i,j-1],MYmax[i,j-1])))
        
        plt.legend(fontsize=9)




        