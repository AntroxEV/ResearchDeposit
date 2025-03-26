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
#fpath3='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_L/DAMP0/'
#fpath4='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Verification/'
fpath5='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/OPENSEESPY/DAMP0/'


cases=['5ms','10ms','15ms','20ms','25ms','50ms']
cases2=['12_5ms','17_5ms','30ms','40ms']
eta=[]
MYmax=[]
MYmaxNL=[]
MYmaxL=[]
penal=[]
DEL_OPT=[]
DEL_NL=[]
DEL_NL2=[]
DEL_L=[]
DELn_OPT=[]
DELn=[]
MYmaxNL2=[]


for i in range(0,len(cases)):
    fname='OptRes_'+cases[i]+'_DEL.obj'
    file_pi = open(fpath + fname, 'rb')
    res=pcle.load(file_pi)
    if  hasattr(res.x, "__len__"):
        tmp=res.x
        eta.append(tmp[0])
    else:
        eta.append(res.x)
    penal.append(res.fun)
    #print(res)
    fname='MYmud_'+cases[i]+'_DEL.txt'  #OPTIMIZED ANALYSES
    AA=np.genfromtxt(fpath + fname)
    MYmax.append(np.max(np.abs(AA[12000:])))
    array_out = rf.rainflow(AA[12000:])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_OPT.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    fname='MYmud_'+cases[i]+'.txt'      #NONLINEAR ANALYSES
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))
    array_out = rf.rainflow(AA[12000:,1])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_NL.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    DELn_OPT.append(np.divide(DEL_OPT[i],DEL_NL[i]))

    fname='MYmud_'+cases[i]+'.txt'     #0 DAMPING ANALYSES
    CC=np.genfromtxt(fpath5 + fname)
    MYmaxL.append(np.max(np.abs(CC[12000:])))  
    array_out = rf.rainflow(CC[12000:])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_L.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    DELn.append(np.divide(DEL_L[i],DEL_NL[i]))

for i in range(0,len(cases2)):
    fname='MYmud_'+cases[i]+'.txt'      #NONLINEAR ANALYSES
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL2.append(np.max(np.abs(AA[12000:,1])))
    array_out = rf.rainflow(AA[12000:,1])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_NL2.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))


    
xs=[5.,10.,15.,20.,25,50]
plt.figure(1)
plt.title('linear fitting')
plt.plot(xs[:-1],eta[:-1],'xk')
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
plt.title('quadratic fitting')
plt.plot(xs,eta,'ok')
xs2=range(0,50)
line = np.multiply(res[2],np.power(xs2,2))+np.multiply(res[1],xs2)+res[0]
plt.plot(xs2, line, '--r', label=r'Fitted line $\eta={:.4f}w_s^2 + {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res[2],res[1],res[0],R_sq ))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel('Added Damping Factor [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,5])


plt.figure(1)
xs3=range(0,25)
slope, intercept, r_value, p_value, std_err = stats.linregress(xs[:-1],eta[:-1])
line = np.multiply(slope,xs3)+intercept
plt.plot(xs3, line, '--r', label=r'Fitted line $\eta={:.2f}w_s {:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))
plt.legend(fontsize=9)

plt.figure(2)
plt.title('bending moment')
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmax,'+r',label='Optimal damping - Linear')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9)

plt.figure(5)
plt.plot(xs,DEL_NL,'ok',label='Nonlinear')
plt.plot(xs,DEL_OPT,'+r',label='Optimal damping - Linear')
plt.plot(xs,DEL_L,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(6)
plt.plot(xs,DELn,'xb',label = 'Normalized DEL - 0-damping ')
plt.plot(xs,DELn_OPT,'+r',label = 'Normalized DEL - Optimal damping ')
plt.plot([xs[0],xs[-1]],[1,1],':k')
#plt.ylim([0,1.6])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(4)

plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmax,'+r',label='Optimal damping - Linear')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(loc='lower right',fontsize=9)




###--------------------------------------------
xs3=[5.,10.,20,50]
eta3=[eta[0], eta[1], eta[3],eta[-1]]
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