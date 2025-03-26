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


fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/Random/'
#fpath3='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_L/DAMP0/'
fpath4='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Verification_Random/'
fpath5='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/OPENSEESPY/DAMP0/Random/'

cases=['_5ms','_10ms','_12_5ms','_15ms','_17_5ms','_20ms','_25ms','_40ms','_50ms']





ncw=len(cases)
eta=np.zeros(shape=(ncw,10))
penal=np.zeros(shape=(ncw,10))
MYmaxV=np.zeros(shape=(ncw,10))
MYmaxVD=np.zeros(shape=(ncw,10))
MYmaxNL=np.zeros(shape=(ncw,10))
MYmaxL=np.zeros(shape=(ncw,10))
DEL_NL=np.zeros(shape=(ncw,10))
DEL_VER=np.zeros(shape=(ncw,10))
DELn=np.zeros(shape=(ncw,10))
DELn_VER=np.zeros(shape=(ncw,10))
DEL_DELVER=np.zeros(shape=(ncw,10))
DELn_DELVER=np.zeros(shape=(ncw,10))
DEL_L=np.zeros(shape=(ncw,10))
DELn_L=np.zeros(shape=(ncw,10))

for i in range(0,len(cases)):
    for j in range(1,11):
     #NONLINEAR ANALYSES 
        fname='MYmud'+cases[i]+'_SET'+str(j)+'.txt'
        AA=np.genfromtxt(fpath2 +cases[i]+'/'+ fname)
        MYmaxNL[i][j-1]=np.max(np.abs(AA[12000:,1]))
        #array_out = rf.rainflow(AA[12000:,1])
        #array_out = array_out[:,array_out[0,:].argsort()]
        #DEL_NL[i][j-1]=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
        rngT=[]
        countT=[]
        for rng, a , count, b,c in rf.extract_cycles(AA[12000:,1]):
            rngT.append(rng)
            countT.append(count)
        DEL_NL[i][j-1]= np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
        

        
        if cases[i] in ('_12_5ms','_17_5ms') :
            fname='MYmud' +'_SET'+str(j)+cases[i]+'.txt'     #0 DAMPING ANALYSES
        else:
            fname='MYmud' +'_SET'+str(j)+cases[i]+'.txt'     #0 DAMPING ANALYSES
        CC=np.genfromtxt(fpath5 + fname)
        MYmaxL[i][j-1]=np.max(np.abs(CC[12000:]))
        #array_out = rf.rainflow(CC[12000:])
        #array_out = array_out[:,array_out[0,:].argsort()]
        #DEL_L[i][j-1]=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
        rngT=[]
        countT=[]
        for rng, a , count, b,c in rf.extract_cycles(CC[12000:]):
            rngT.append(rng)
            countT.append(count)
        DEL_L[i][j-1]= np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
        DELn_L[i][j-1]=np.divide(DEL_L[i][j-1],DEL_NL[i][j-1])

        fname='MYmud'+'_SET'+str(j)+cases[i]+'_OPTVER2.txt'
        AA=np.genfromtxt(fpath4 + fname)
        MYmaxV[i][j-1]=np.max(np.abs(AA[12000:]))
        #array_out = rf.rainflow(AA[12000:])
        #array_out = array_out[:,array_out[0,:].argsort()]
        #DEL_VER[i][j-1]=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
        rngT=[]
        countT=[]
        for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
            rngT.append(rng)
            countT.append(count)
        DEL_VER[i][j-1]= np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
        DELn_VER[i][j-1]=(np.divide(DEL_VER[i][j-1],DEL_NL[i][j-1]))
    
        fname='MYmud'+'_SET'+str(j)+cases[i]+'_DELVER2.txt'
        AA=np.genfromtxt(fpath4 + fname)
        MYmaxVD[i][j-1]=np.max(np.abs(AA[12000:]))
        #array_out = rf.rainflow(AA[12000:])
        #array_out = array_out[:,array_out[0,:].argsort()]
        #DEL_DELVER[i][j-1]=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
        rngT=[]
        countT=[]
        for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
            rngT.append(rng)
            countT.append(count)
        DEL_DELVER[i][j-1]= np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
        DELn_DELVER[i][j-1]=(np.divide(DEL_DELVER[i][j-1],DEL_NL[i][j-1]))
    
xs=[5.,10.,12.5,15.,17.5,20.,25,40,50]



plt.figure(1,figsize=(5, 4))
plt.plot(xs,MYmaxNL,'ok',label='__nolegend__')
plt.plot(xs,MYmaxV,'+r',label='__nolegend__')
plt.plot(xs,MYmaxVD,'.c',label='__nolegend__')
plt.plot(xs,MYmaxL,'xb',label='__nolegend__')
plt.plot(xs[-1],MYmaxNL[-1][0],'ok',mfc='none',label='Nonlinear Model')
plt.plot(xs[-1],MYmaxV[-1][0],'+r',label='costant damping - Max Bending')
plt.plot(xs[-1],MYmaxVD[-1][0],'.c',label='costant damping - DEL')
plt.plot(xs[-1],MYmaxL[-1][0],'xb',label='0-damping')
plt.plot(xs,np.mean(MYmaxV,1),'--r')
plt.plot(xs,np.mean(MYmaxVD,1),'--c')
plt.plot(xs,np.mean(MYmaxL,1),'--b')
plt.plot(xs,np.mean(MYmaxNL,1),':k')
#plt.plot(xs,MYmaxL,'xb',label='0-damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9,loc='lower right')

plt.figure(2,figsize=(5, 4))
plt.plot(xs,DEL_NL,'ok',label='__nolegend__')
plt.plot(xs,DEL_VER,'+r',label='__nolegend__')
plt.plot(xs,DEL_DELVER,'.c',label='__nolegend__')
plt.plot(xs,DEL_L,'xb',label='__nolegend__')
plt.plot(xs[-1],DEL_NL[-1][0],'ok',mfc='none',label='Nonlinear Model')
plt.plot(xs[-1],DEL_VER[-1][0],'+r',label='constant damping - Max Bending')
plt.plot(xs[-1],DEL_DELVER[-1][0],'.c',label='costant damping - DEL')
plt.plot(xs[-1],DEL_L[-1][0],'xb',label='0-damping')
plt.plot(xs,np.mean(DEL_VER,1),'--r')
plt.plot(xs,np.mean(DEL_DELVER,1),'--c')
plt.plot(xs,np.mean(DEL_L,1),'--b')
plt.plot(xs,np.mean(DEL_NL,1),':k')
#plt.plot(xs,DEL_L,'xb',label='0-damping - Linear')
#plt.plot(xs,DEL_DELVER,'.c',label='Fitted damping - DEL')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [kNm]',**hfont, fontsize=12)
#plt.legend(fontsize=9)


plt.figure(3,figsize=(5, 4))
#plt.plot(xs,DELn,'xb',label = '0-damping')
plt.plot(xs,DELn_VER,'+r',label = 'constant damping - Max Bending ')
plt.plot(xs,DELn_DELVER,'.c',label = 'constant damping - DEL ')
plt.plot(xs,DELn_L,'xb',label = '0-damping ')
plt.plot([xs[0],xs[-1]],[1,1],':k')
plt.plot(xs,np.mean(DELn_VER,1),'--r')
plt.plot(xs,np.mean(DELn_L,1),'--b')
plt.plot(xs,np.mean(DELn_DELVER,1),'--c')
plt.ylim([0,2.0])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
#plt.legend(fontsize=9)


plt.figure(10,figsize=(5, 4))
rat=np.divide(np.mean(MYmaxV,1),np.mean(MYmaxNL,1))
plt.plot(xs,rat,'--k',marker='+',label = ' Bending Moment - OPT:Max Bending ')
rat=np.divide(np.mean(MYmaxVD,1),np.mean(MYmaxNL,1))
plt.plot(xs,rat,'--c',marker='+',label = ' Bending Moment - OPT:DEL ')
plt.plot([xs[0],xs[-1]],[1,1],':k')
#rat=np.divide(np.mean(DEL_VER,1),np.mean(DEL_NL,1))
#plt.plot(xs,rat,'--k',marker='.',label = 'DEL - OPT:Max Bending ')
#rat=np.divide(np.mean(DEL_DELVER,1),np.mean(DEL_NL,1))
#plt.plot(xs,rat,'--c',marker='.',label = 'DEL - OPT:Max Bending ')
#hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized Peak Bending Moment [-]',**hfont, fontsize=12)
plt.locator_params(axis='both', nbins=15)
plt.ylim([.5,1.5])
plt.legend(fontsize=9)

plt.figure(11,figsize=(5, 4))
#rat=np.divide(np.mean(MYmaxV,1),np.mean(MYmaxNL,1))
#plt.plot(xs,rat,'--k',marker='+',label = ' Bending Moment - OPT:Max Bending ')
#rat=np.divide(np.mean(MYmaxVD,1),np.mean(MYmaxNL,1))
#plt.plot(xs,rat,'--c',marker='+',label = ' Bending Moment - OPT:DEL ')
rat=np.divide(np.mean(DEL_VER,1),np.mean(DEL_NL,1))
plt.plot(xs,rat,'--k',marker='.',label = 'DEL - OPT:Max Bending ')
plt.plot([xs[0],xs[-1]],[1,1],':k')
rat=np.divide(np.mean(DEL_DELVER,1),np.mean(DEL_NL,1))
plt.plot(xs,rat,'--c',marker='.',label = 'DEL - OPT:DEL ')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.locator_params(axis='both', nbins=15)
plt.ylim([.5,1.5])
plt.legend(fontsize=9)