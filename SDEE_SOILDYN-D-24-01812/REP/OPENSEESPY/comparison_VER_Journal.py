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
MYmaxV2=[]
MYmaxVD=[]
MYmaxVD2=[]
DEL_VER=[]
DELn_VER=[]
DEL_VER2=[]
DELn_VER2=[]

DEL_DELVER=[]
DELn_DELVER=[]
DEL_DELVER2=[]
DELn_DELVER2=[]
    
    
for i in range(0,len(cases)):
    fname='MYmud_'+cases[i]+'.txt'      #NONLINEAR ANALYSES
    AA=np.genfromtxt(fpath2 + fname)
    MYmaxNL.append(np.max(np.abs(AA[12000:,1])))
    #array_out = rf.rainflow(AA[12000:,1])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_NL.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:,1]):
        rngT.append(rng)
        countT.append(count)
    DEL_NL.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))

    fname='MYmud_'+cases[i]+'.txt'     #0 DAMPING ANALYSES
    CC=np.genfromtxt(fpath5 + fname)
    MYmaxL.append(np.max(np.abs(CC[12000:])))  
    #array_out = rf.rainflow(CC[12000:])
    #rray_out = array_out[:,array_out[0,:].argsort()]
    #DEL_L.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(CC[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_L.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn.append(np.divide(DEL_L[i],DEL_NL[i]))

    fname='MYmud_'+cases[i]+'_OPTVER.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxV.append(np.max(np.abs(AA[12000:])))
    #array_out = rf.rainflow(AA[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_VER.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_VER.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_VER.append(np.divide(DEL_VER[i],DEL_NL[i]))
    
    fname='MYmud_'+cases[i]+'_DELVER.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxVD.append(np.max(np.abs(AA[12000:])))
    #array_out = rf.rainflow(AA[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_DELVER.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_DELVER.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_DELVER.append(np.divide(DEL_DELVER[i],DEL_NL[i]))
    
    
    fname='MYmud_'+cases[i]+'_OPTVER2.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxV2.append(np.max(np.abs(AA[12000:])))
    #array_out = rf.rainflow(AA[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_VER.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_VER2.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_VER2.append(np.divide(DEL_VER2[i],DEL_NL[i]))
    
    fname='MYmud_'+cases[i]+'_DELVER2.txt'
    AA=np.genfromtxt(fpath4 + fname)
    MYmaxVD2.append(np.max(np.abs(AA[12000:])))
    #array_out = rf.rainflow(AA[12000:])
    #array_out = array_out[:,array_out[0,:].argsort()]
    #DEL_DELVER.append(np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef))
    rngT=[]
    countT=[]
    for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
        rngT.append(rng)
        countT.append(count)
    DEL_DELVER2.append(np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef))
    DELn_DELVER2.append(np.divide(DEL_DELVER2[i],DEL_NL[i]))
    
xs=[5.,10.,12.5,15.,17.5,20.,25,30,40,50]



plt.figure(1,figsize=(5, 4))
plt.plot(xs,MYmaxNL,'ok',label='Nonlinear')
plt.plot(xs,MYmaxV,'+r',label='constant damping - Max Moment')
plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
plt.plot(xs,MYmaxVD,'.c',label='constant damping - DEL')
plt.plot(xs,MYmaxV2,'1m',label='linear damping - Max Moment')
plt.plot(xs,MYmaxVD2,'2g', mfc='none',label='linear damping - DEL')

hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
#plt.legend(fontsize=9,loc='lower right')

plt.figure(2,figsize=(5, 4))
plt.plot(xs,DEL_NL,'ok',label='Nonlinear')
plt.plot(xs,DEL_VER,'+r',label='constant damping - Max Moment')
plt.plot(xs,DEL_L,'xb',label='0-damping - Linear')
plt.plot(xs,DEL_DELVER,'.c',label='constant damping - DEL')
plt.plot(xs,DEL_VER2,'1r',label='linear damping - Max Moment')
plt.plot(xs,DEL_DELVER2,'oc', mfc='none',label='linear damping - DEL')

hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(3,figsize=(5, 4))
plt.plot(xs,DELn,'xb',label = '0-damping')
plt.plot(xs,DELn_VER,'+r',label = 'constant damping - Max Moment ')
plt.plot(xs,DELn_DELVER,'.c',label = 'constant damping - DEL ')
plt.plot(xs,DELn_VER2,'1m',label='linear damping - Max Moment')
plt.plot(xs,DELn_DELVER2,'2g', mfc='none',label='linear damping - DEL')
plt.plot([xs[0],xs[-1]],[1,1],':k')
plt.ylim([0,2.5])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


