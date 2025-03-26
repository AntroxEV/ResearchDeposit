# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 08:15:55 2023

@author: at924
"""

import pickle as pcle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

fpath='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR_OPT2/'
fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/'
cases=['5ms','12_5ms','15ms','20ms','50ms']
eta1=[]
eta2=[]
MYmax=[]
MYmaxNL=[]
penal=[]
for i in range(0,len(cases)):
    fname='OptRes2_'+cases[i]+'.obj'
    file_pi = open(fpath + fname, 'rb')
    res=pcle.load(file_pi)
    eta1.append(res.x[0])
    eta2.append(res.x[1])
    penal.append(res.fun)
    #print(res)
    fname='MYmud_'+cases[i]+'_OPT2.txt'
    AA=np.genfromtxt(fpath + fname)
    MYmax.append(np.max(np.abs(AA[12000:])))
    fname='MYmud_'+cases[i]+'.txt'
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))


plt.figure(1)
xs=[5.,12.5,15.,20.,50]
plt.plot(xs,eta1,'*k')
plt.plot(xs,eta2,'+r')
#.xticks(xs,cases)

A=np.zeros([len(xs),2])
A[:,0]=[0 for _ in range(0,len(xs))]
A[:,1]=np.transpose(xs)

res=np.dot(np.linalg.pinv(A),eta1)
line = np.multiply(res[1],xs)+res[0]

corr_matrix = np.corrcoef(eta1, line)
corr = corr_matrix[0,1]
R_sq = corr**2

plt.figure(3)
xs=[5.,12.5,15.,20.,50]
plt.plot(xs,eta1,'*')
plt.plot(xs, line, '--r', label=r'Fitted line $\eta={:.2f}w_s {:.2f}$,  $R^2 =$ {:.2f}'.format(res[1],res[0],R_sq ))
plt.legend(fontsize=9)


plt.figure(1)
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,eta1)
line = np.multiply(slope,xs)+intercept



plt.plot(xs, line, '--r', label=r'Fitted line $\eta={:.2f}w_s {:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))

plt.legend(fontsize=9)

plt.figure(2)
plt.plot(xs,MYmax,'ok')
plt.plot(xs,MYmaxNL,'+r')

