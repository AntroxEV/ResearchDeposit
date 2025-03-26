# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 08:15:55 2023

@author: at924
"""

import pickle as pcle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import rainflow as rf

Nref=2E8;
mcoef=4;

fpath='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/'
fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/'
fpath3='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/OPENSEESPY/DAMP0/'


cases=['5ms','10ms','15ms','20ms','25ms','30ms','40ms','50ms']
#cases3=['5ms','10ms','15ms','20ms','25ms','50ms']

eta=[]
MYmax=[]
MYmaxNL=[]
MYmaxL=[]
penal=[]
DEL_OPT=[]
DEL_NL=[]
DEL_L=[]
DELn_OPT=[]
DELn_L=[]
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
    #array_out = rf.rainflow(AA[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_OPT.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    #DEL_OPT.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_OPT.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    
    fname='MYmud_'+cases[i]+'.txt'
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))
    #array_out = rf.rainflow(AA[12000:,1])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_NL.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    #DEL_NL.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:,1]):
        rngT.append(rng)
        countT.append(count)
    DEL_NL.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_OPT.append(np.divide(DEL_OPT[i],DEL_NL[i]))

    fname='MYmud_'+cases[i]+'.txt'
    CC=np.genfromtxt(fpath3 + fname)
    MYmaxL.append(np.max(np.abs(CC[12000:])))
    #array_out = rf.rainflow(CC[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_L.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    #DEL_L.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(CC[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_L.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_L.append(np.divide(DEL_L[i],DEL_NL[i]))
    
    

    
    

plt.figure(1,figsize=(5, 4))
xs=[5.,10.,15.,20.,25,30,40,50]
plt.plot(xs[:-3],eta[:-3],'ok',label='Optimal damping')
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

plt.figure(3,figsize=(5, 4))
plt.plot(xs,eta,'ok',label='Optimal Values')
xs2=range(0,51)
#line = np.multiply(res[2],np.power(xs2,2))+np.multiply(res[1],xs2)+res[0]
#plt.plot(xs2, line, '--r', label=r'Fitted line $\eta={:.4f}w_s^2 + {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res[2],res[1],res[0],R_sq ))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel(r'Stiffness proportional coefficient $\beta$  [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,1])


###--------------------------------------------
xs3=[10.,20,40]
eta3=[ eta[1], eta[3],eta[-2]]
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
plt.plot(50,0.2,'ob',label='$2^{nd}$ local optimum')
xs2=range(0,51)
#line = np.multiply(res3[2],np.power(xs2,2))+np.multiply(res3[1],xs2)+res3[0]
#plt.plot(xs2, line, '--r', label=r' $\eta={:.4f}w_s^2 {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res3[2],res3[1],res3[0],R_sq ))
plt.plot([0,50],[np.mean(eta[:-1]),np.mean(eta[:-1])],'--r',label=r' $\beta_{mean}$' +'={:.4f}'.format(np.mean(eta[:-1])))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel(r'Stiffness proportional coefficient $\beta$  [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,5])



plt.figure(2,figsize=(5, 4))
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear Model')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear Model')
plt.plot(xs,MYmax,'+r',label='Optimal Damping - Linear Model')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(loc='lower right',fontsize=9)



plt.figure(1,figsize=(5, 4))
slope, intercept, r_value, p_value, std_err = stats.linregress(xs[:-3],eta[:-3])
line = np.multiply(slope,xs[:-3])+intercept
plt.plot(xs[:-3], line, '--r', label=r' $\beta={:.2f}w_s {:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel(r'Stiffness proportional coefficient $\beta$  [-]',**hfont, fontsize=12)
#plt.xlim([0,25])
plt.legend(fontsize=9)


plt.figure(7,figsize=(5, 4))
plt.plot(xs,DELn_L,'xb',label = ' 0-damping - Linear Model ')
plt.plot(xs,DELn_OPT,'+r',label = 'Optimal damping - Linear Model ')
plt.plot([xs[0],xs[-1]],[1,1],':k',label='Nonlinear Model')
plt.ylim([0,2.5])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
