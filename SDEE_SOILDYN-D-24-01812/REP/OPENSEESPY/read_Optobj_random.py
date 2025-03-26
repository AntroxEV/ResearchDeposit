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


fpath='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Random/'
fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/Random/'
cases=['_5ms','_10ms','_15ms','_20ms','_25ms']
ncw=len(cases)
DEL_NL=np.zeros(shape=(ncw,10))
DEL_OPT=np.zeros(shape=(ncw,10))
DELn=np.zeros(shape=(ncw,10))
eta=np.zeros(shape=(ncw,10))
penal=np.zeros(shape=(ncw,10))
MYmax=np.zeros(shape=(ncw,10))
MYmaxNL=np.zeros(shape=(ncw,10))
for i in range(0,len(cases)):
    for j in range(1,11):
        fname='OptRes'+cases[i]+'_SET'+str(j)+'.obj'
        file_pi = open(fpath + cases[i]+'/'+ fname, 'rb')
        res=pcle.load(file_pi)
        penal[i][j-1]=res.fun
        if res.fun<1e10:
            eta[i][j-1]=res.x
        else:
            eta[i][j-1]=None
        #print(res)
        fname='MYmud_SET'+str(j)+cases[i]+'_OPT.txt'
        AA=np.genfromtxt(fpath +cases[i]+'/'+ fname)
        array_out = rf.rainflow(AA[12000:])
        array_out = array_out[:,array_out[0,:].argsort()]
        DEL_OPT[i][j-1]=np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef)
        MYmax[i][j-1]=np.max(np.abs(AA[12000:]))
        fname='MYmud'+cases[i]+'_SET'+str(j)+'.txt'
        AA=np.genfromtxt(fpath2 +cases[i]+'/'+ fname)
        MYmaxNL[i][j-1]=np.max(np.abs(AA[12000:,1]))
        array_out = rf.rainflow(AA[12000:,1])
        array_out = array_out[:,array_out[0,:].argsort()]
        DEL_NL[i][j-1]=np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef)
DELn=(np.divide(DEL_OPT,DEL_NL))
#penal[penal>0.01]=0


plt.figure(1)
xs=[5.,10.,15.,17.5,20.,25.]
plt.plot([xs[0]]*len(eta[0]),eta[0],'*',color='black')
plt.plot([xs[1]]*len(eta[1]),eta[1],'*',color='red')
plt.plot([xs[2]]*len(eta[2]),eta[2],'*',color='blue')
plt.plot([xs[3]]*len(eta[3]),eta[3],'*',color='green')
plt.plot([xs[4]]*len(eta[4]),eta[4],'*',color='magenta')
plt.plot([xs[5]]*len(eta[5]),eta[5],'+',color='magenta')

plt.figure(2)


eta0=np.nan_to_num(eta, copy=True, nan=0.0, posinf=None, neginf=None)
etam=np.mean(eta0,axis=1)
plt.plot(xs,etam,'o',color='black')
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,etam)
line = np.multiply(slope,xs)+intercept
plt.plot(xs, line, '--r', label=r'Fitted line $\eta={:.2f}w_s {:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))

plt.legend(fontsize=9)
# #.xticks(xs,cases)


plt.figure(3)
xs=[5.,10.,15.,17.5,20.,25.]
plt.plot([xs[0]]*len(penal[0]),penal[0],'*',color='black')
plt.plot([xs[1]]*len(penal[1]),penal[1],'*',color='red')
plt.plot([xs[2]]*len(penal[2]),penal[2],'*',color='blue')
plt.plot([xs[3]]*len(penal[3]),penal[3],'*',color='green')
plt.plot([xs[4]]*len(penal[4]),penal[4],'*',color='magenta')
plt.plot([xs[5]]*len(penal[5]),penal[5],'+',color='magenta')
#plt.plot(penal,eta,'*')
# plt.plot(xs,MYmax,'^')

plt.figure(4)
plt.plot(xs,MYmaxNL,'ok',label='__nolegend__')
plt.plot(xs,MYmax,'+r',label='__nolegend__')
plt.plot(xs[-1],MYmaxNL[-1][0],'ok',label='Nonlinear')
plt.plot(xs[-1],MYmax[-1][0],'+r',label='Optimal damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9)
#plt.figure(3)
#plt.plot([xs[i] for i in range(0,4)],[eta[i] for i in range(0,4)])

#plt.figure(4)
#ax1=plt.hist(np.transpose(eta), 6, density=True, histtype='bar', stacked=True)
#ax1.set_title('stacked bar')

plt.figure(5)
plt.hist(eta[0],4,color='k',label='5ms')
plt.hist(eta[1],4,color='b',label='10ms')
plt.hist(eta[2],4,color='r',label='15ms')
plt.hist(eta[3],4,color='m',label='17.5ms')
plt.hist(eta[4],4,color='g',label='20ms')
plt.hist(eta[5],4,color='y',label='25ms')
plt.legend()



plt.figure(6)
plt.plot(xs,DEL_NL,'ok',label='__nolegend__')
plt.plot(xs,DEL_OPT,'+r',label='__nolegend__')
plt.plot(xs[-1],DEL_NL[-1][0],'ok',label='Nonlinear')
plt.plot(xs[-1],DEL_OPT[-1][0],'+r',label='Optimal damping - Linear')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(1001)
[plt.plot(np.repeat(xs[i], 10),DELn[i],'+r',label = '__nolegend__') for i in range(0,ncw)]
plt.plot([xs[0],xs[-1]],[1,1],':k')
#plt.ylim([0,1.6])
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.legend()
