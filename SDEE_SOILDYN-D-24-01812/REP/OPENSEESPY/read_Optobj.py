# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 08:15:55 2023

@author: at924
"""

import pickle as pcle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

fpath='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/'
fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/'
fpath3='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_L/DAMP0/'
fpath4='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Verification/'
fpath5='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/OPENSEESPY/DAMP0/'

cases=['5ms','10ms','15ms','20ms','25ms','50ms']
cases3=['5ms','10ms','15ms','20ms','25ms','50ms']
caseV=['30ms','40ms']
eta=[]
MYmax=[]
MYmaxNL=[]
MYmaxL=[]
penal=[]
MYmaxV=[]
MYmaxNLV=[]
MYmaxLV=[]
for i in range(0,len(cases)):
    fname='OptRes_'+cases[i]+'.obj'
    file_pi = open(fpath + fname, 'rb')
    res=pcle.load(file_pi)
    if  hasattr(res.x, "__len__"):
        tmp=res.x
        eta.append(tmp[0])
    else:
        eta.append(res.x)
    penal.append(res.fun)
    #print(res)
    fname='MYmud_'+cases[i]+'_OPT.txt'
    AA=np.genfromtxt(fpath + fname)
    MYmax.append(np.max(np.abs(AA[12000:])))
    fname='MYmud_'+cases[i]+'.txt'
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))

for i in range(0,len(cases)):
    fname='MYmud_'+cases[i]+'.txt'
    CC=np.genfromtxt(fpath3 + fname)
    MYmaxL.append(np.max(np.abs(CC[12000:,1])))  
    
    
for i in range(0,len(caseV)):
    fname='MYmud_'+caseV[i]+'.txt'
    CC=np.genfromtxt(fpath2 + fname)
    MYmaxNLV.append(np.max(np.abs(CC[12000:,1])))
    fname='MYmud_'+caseV[i]+'_VER.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxV.append(np.max(np.abs(AA[12000:])))
    fname='MYmud_'+caseV[i]+'.txt'
    AA=np.genfromtxt(fpath5 + fname)
    MYmaxLV.append(np.max(np.abs(AA[12000:])))
    
    

plt.figure(1)
xs=[5.,10.,12.5,15.,17.5,20.,25,50]
plt.plot(xs[:-1],eta[:-1],'ok',label='Optimal damping')
#.xticks(xs,cases)

A=np.zeros([len(xs),3])
A[:,0]=[1 for _ in range(0,len(xs))]
A[:,1]=np.transpose(xs)
A[:,2]=np.power(np.transpose(xs),2)

res=np.dot(np.linalg.pinv(A),eta)
line = np.multiply(res[2],np.power(xs,2))+np.multiply(res[1],xs)+res[0]

corr_matrix = np.corrcoef(eta, line)
corr = corr_matrix[0,1]
R_sq = corr**2

plt.figure(3)
xs=[5.,10.,12.5,15.,17.5,20.,25,50]
plt.plot(xs,eta,'ok')
xs2=range(0,51)
line = np.multiply(res[2],np.power(xs2,2))+np.multiply(res[1],xs2)+res[0]
plt.plot(xs2, line, '--r', label=r'Fitted line $\eta={:.4f}w_s^2 + {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res[2],res[1],res[0],R_sq ))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel('Added Damping Factor [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,5])


###--------------------------------------------
xs3=[5.,10.,20,50]
eta3=[eta[0], eta[1], eta[5],eta[7]]
A=np.zeros([len(xs3),3])
A[:,0]=[1 for _ in range(0,len(xs3))]
A[:,1]=np.transpose(xs3)
A[:,2]=np.power(np.transpose(xs3),2)

res3=np.dot(np.linalg.pinv(A),eta3)
line3 = np.multiply(res3[2],np.power(xs3,2))+np.multiply(res3[1],xs3)+res3[0]

corr_matrix = np.corrcoef(eta3, line3)
corr = corr_matrix[0,1]
R_sq = corr**2


plt.figure(3)
plt.plot(xs3,eta3,'ob')
xs2=range(0,51)
line = np.multiply(res3[2],np.power(xs2,2))+np.multiply(res3[1],xs2)+res3[0]
plt.plot(xs2, line, '--k', label=r'Fitted line $\eta={:.4f}w_s^2 + {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res3[2],res3[1],res3[0],R_sq ))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel('Added Damping Factor [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,5])



plt.figure(2)
xs=[5.,10.,12.5,15.,17.5,20.,25,50]
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(loc='lower right',fontsize=9)

plt.figure(4)
xs=[5.,10.,12.5,15.,17.5,20.,25,50]
xs4=[30,40]
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmax,'+r',label='Optimal damping - Linear')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
plt.plot(xs4,MYmaxNLV,'or',label='Nonlinear')
plt.plot(xs4,MYmaxV,'+r')
plt.plot(xs4,MYmaxLV,'xb')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(loc='lower right',fontsize=9)

plt.figure(1)
slope, intercept, r_value, p_value, std_err = stats.linregress(xs[:-1],eta[:-1])
line = np.multiply(slope,xs[:-1])+intercept
plt.plot(xs[:-1], line, '--k', label=r'Fitted line $\eta={:.2f}w_s {:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel('Added Damping Factor [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


