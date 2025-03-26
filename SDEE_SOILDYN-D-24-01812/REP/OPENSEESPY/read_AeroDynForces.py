# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:47:46 2024

@author: at924
"""

import numpy as np
import matplotlib.pyplot as plt


#fpath='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/AerodynForces/Random/'
cases=['_5ms','_10ms','_15ms','_20ms','_25ms']

cut=24000
time=np.linspace(cut*0.005, 720,144001-cut)
# FXmax=np.zeros(shape=(5,10))

# for i in range(0,len(cases)):
#     for j in range(1,11):
#         data=np.genfromtxt(fpath + cases[i] + '/SET' + str(j) + '/FX' +cases[i] + '.txt')
#         plt.figure(i)
#         plt.plot(time,data[cut:,1],label=str(j))
#         FXmax[i,j-1]=np.max(np.abs(data[cut:,1]))
        

# plt.figure(100)
# xs=[5.,10.,15.,20.,25.]
# plt.plot([xs[0]]*10,FXmax[0,:],'*',color='black')
# plt.plot([xs[1]]*10,FXmax[1,:],'*',color='red')
# plt.plot([xs[2]]*10,FXmax[2,:],'*',color='blue')
# plt.plot([xs[3]]*10,FXmax[3,:],'*',color='green')
# plt.plot([xs[4]]*10,FXmax[4,:],'*',color='magenta')




fpath='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/AerodynForces/'
ls=['5m/s','15m/s','25m/s']
cl=['k','--r',':b']
ind=0
for i in range(0,len(cases),2):
        data=np.genfromtxt(fpath + '/FX' +cases[i] + '.txt')
        plt.figure(200,figsize=(5, 5))
        plt.plot(time-120,data[cut:,1]/1000,cl[ind],label=ls[ind])
        ind+=1

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Fore-aft force Fy [kN]',fontsize=14)
plt.xlabel('time [s]',fontsize=16)
plt.legend(fontsize=12)

ind=0
for i in range(0,len(cases),2):
        data=np.genfromtxt(fpath + '/MY' +cases[i] + '.txt')
        plt.figure(300,figsize=(5, 5))
        plt.plot(time-120,data[cut:,1]/1000,cl[ind],label=ls[ind])
        ind+=1

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Moment Mr [kNm]',fontsize=14)
plt.xlabel('time [s]',fontsize=16)
plt.legend(fontsize=12)


ind=0
for i in range(0,len(cases),2):
        data=np.genfromtxt(fpath + '/FZ' +cases[i] + '.txt')
        plt.figure(400,figsize=(5, 5))
        plt.plot(time-120,data[cut:,1]/1000,cl[ind],label=ls[ind])
        ind+=1

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Vertical force Fz [kNm]',fontsize=14)
plt.xlabel('time [s]',fontsize=16)
plt.legend(fontsize=12)