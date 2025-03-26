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


fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/'
#fpath3='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_L/DAMP0/'
fpath4='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Verification/'
fpath5='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/OPENSEESPY/DAMP0/'

cases=['5ms','10ms','12_5ms','15ms','17_5ms','20ms','25ms','30ms','40ms','50ms']

nC=len(cases)
eta=[]

MYmaxNL=[]
MYmaxL=[]
penal=[]
DEL_NL=[]
DEL_L=[]
DELn_OPT=[]
DELn=[]
MYmaxV=[]
DEL_VER=[]
DELn_VER=[]



    
    
for i in range(0,len(cases)):
    fname='MYmud_'+cases[i]+'.txt'      #NONLINEAR ANALYSES
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))
    array_out = rf.rainflow(AA[12000:,1])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_NL.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))

    fname='MYmud_'+cases[i]+'.txt'     #0 DAMPING ANALYSES
    CC=np.genfromtxt(fpath5 + fname)
    MYmaxL.append(np.max(np.abs(CC[12000:])))  
    array_out = rf.rainflow(CC[12000:])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_L.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    DELn.append(np.divide(DEL_L[i],DEL_NL[i]))

    fname='MYmud_'+cases[i]+'_OPTVER.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxV.append(np.max(np.abs(AA[12000:])))
    array_out = rf.rainflow(AA[12000:])
    array_out = array_out[:,array_out[0,:].argsort()]
    DEL_VER.append(np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef))
    DELn_VER.append(np.divide(DEL_VER[i],DEL_NL[i]))
    
xs=[5.,10.,12.5,15.,17.5,20.,25,30,40,50]



plt.figure(1)
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmaxV,'+r',label='Fitted damping - Linear')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9)

plt.figure(2)
plt.plot(xs,DEL_NL,'ok',label='Nonlinear')
plt.plot(xs,DEL_VER,'+r',label='Fitted damping - Linear')
plt.plot(xs,DEL_L,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(3)
plt.plot(xs,DELn,'xb',label = 'Normalized DEL - 0-damping ')
plt.plot(xs,DELn_VER,'+r',label = 'Normalized DEL - Fitted damping ')
plt.plot([xs[0],xs[-1]],[1,1],':k')
#plt.ylim([0,1.6])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


