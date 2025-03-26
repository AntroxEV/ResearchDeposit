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
xs=[5.,10.,12.5,15.,17.5,20.,25.,40,50]
fpath='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Random_DEL/'
fpath2='C:/Users/at924\OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Results_Stage1/SEISMO_NL/Random/'
fpath3='C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/Random/'
cases=['_5ms','_10ms','_12_5ms','_15ms','_17_5ms','_20ms','_25ms','_40ms','_50ms']

ncw=len(cases)
eta=np.zeros(shape=(ncw,10))
penal=np.zeros(shape=(ncw,10))
MYmax=np.zeros(shape=(ncw,10))
MYmaxNL=np.zeros(shape=(ncw,10))
DEL_NL=np.zeros(shape=(ncw,10))
DEL_OPT=np.zeros(shape=(ncw,10))
DELn=np.zeros(shape=(ncw,10))


for i in range(0,len(cases)):
    for j in range(1,11):
        fname='OptRes'+cases[i]+'_SET'+str(j)+'_DEL.obj'
        file_pi = open(fpath + cases[i]+'/'+ fname, 'rb')
        res=pcle.load(file_pi)
        penal[i][j-1]=res.fun
        eta[i][j-1]=res.x
        fname='MYmud_SET'+str(j)+cases[i]+'_DEL.txt'
        AA=np.genfromtxt(fpath +cases[i]+'/'+ fname)
        MYmax[i][j-1]=np.max(np.abs(AA[12000:]))
        rngT=[]
        countT=[]
        for rng, a , count, b,c in rf.extract_cycles(AA[12000:]):
            rngT.append(rng)
            countT.append(count)
        DEL_OPT[i][j-1] = np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
        #array_out = rf.rainflow(AA[12000:])
        #array_out = array_out[:,array_out[0,:].argsort()]
        #DEL_OPT[i][j-1]=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
        
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
        DELn[i][j-1]=(np.divide(DEL_OPT[i][j-1],DEL_NL[i][j-1]))


A=np.zeros([len(xs),3])
A[:,0]=[1 for _ in range(0,len(xs))]
A[:,1]=np.transpose(xs)
A[:,2]=np.power(np.transpose(xs),2)

res3=np.dot(np.linalg.pinv(A),np.mean(eta,1))
line = np.multiply(res3[2],np.power(xs,2))+np.multiply(res3[1],xs)+res3[0]

corr_matrix = np.corrcoef(np.mean(eta,1), line)
corr = corr_matrix[0,1]
R_sq = corr**2


xs=[5.,10.,12.5,15.,17.5,20.,25.,40,50]
plt.figure(1,figsize=(5, 4))
plt.plot(xs,eta,'ok',mfc='none',label='__nolegend__')
plt.plot(xs[-1],eta[0][-1],'ok',mfc='none',label='Optimal damping - DEL')
#plt.plot(xs,np.mean(eta,1),'--k',label='__nolegend__')
xs2=range(0,51)
#line = np.multiply(res3[2],np.power(xs2,2))+np.multiply(res3[1],xs2)+res3[0]
#plt.plot(xs2, line, '--r', label=r' $\eta={:.4f}w_s^2 +{:.4f}w_s  {:.4f}$,  $R^2 =$ {:.2f}'.format(res3[2],res3[1],res3[0],R_sq ))
plt.plot(xs, [np.mean(eta)]*len(xs), '--k', label=r'$\beta_{mean}$' + '={:.4f}'.format(np.mean(eta)))
slope, intercept, r_value, p_value, std_err = stats.linregress(xs,np.mean(eta,1))
line = np.multiply(slope,xs)+intercept
plt.plot(xs, line, '--r', label=r' $\beta={:.2f}w_s +{:.2f}$, $R^2 =$ {:.2f} '.format(slope,intercept,r_value))
#line = np.multiply(res[2],np.power(xs2,2))+np.multiply(res[1],xs2)+res[0]
#plt.plot(xs2, line, '--r', label=r'Fitted line $\eta={:.4f}w_s^2 + {:.4f}w_s + {:.4f}$,  $R^2 =$ {:.2f}'.format(res[2],res[1],res[0],R_sq ))
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ylabel(r'Stiffness proportional coefficient $\beta$  [-]',**hfont, fontsize=12)
plt.legend(fontsize=9)
plt.ylim([0,1])







plt.figure(3)
plt.title('penalty')
plt.plot([xs[0]]*len(penal[0]),penal[0],'*',color='black')
plt.plot([xs[1]]*len(penal[1]),penal[1],'*',color='red')
plt.plot([xs[2]]*len(penal[2]),penal[2],'*',color='blue')
plt.plot([xs[3]]*len(penal[3]),penal[3],'*',color='green')
plt.plot([xs[4]]*len(penal[4]),penal[4],'*',color='magenta')
#plt.plot(penal,eta,'*')
# plt.plot(xs,MYmax,'^')

plt.figure(4,figsize=(5, 4))
plt.plot(xs,MYmaxNL,'ok',mfc='none',label='__nolegend__')
plt.plot(xs,MYmax,'+r',label='__nolegend__')
plt.plot(xs[-1],MYmaxNL[-1][0],'ok',mfc='none',label='Nonlinear Model')
plt.plot(xs[-1],MYmax[-1][0],'+r',label='Optimal damping - Linear')
plt.plot(xs,np.mean(MYmax,1),'--r')
plt.plot(xs,np.mean(MYmaxNL,1),':k')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('absolute maximum bending moment [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9)
#plt.figure(3)
#plt.plot([xs[i] for i in range(0,4)],[eta[i] for i in range(0,4)])



plt.figure(6)
plt.plot(xs,DEL_NL,'ok',mfc='none',label='__nolegend__')
plt.plot(xs,DEL_OPT,'+r',label='__nolegend__')
plt.plot(xs[-1],DEL_NL[-1][0],'ok',mfc='none',label='Nonlinear')
plt.plot(xs[-1],DEL_OPT[-1][0],'+r',label='Optimal damping - Linear')
plt.plot(xs,np.mean(DEL_OPT,1),'--r')
plt.plot(xs,np.mean(DEL_NL,1),':k')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('DEL [kNm]',**hfont, fontsize=12)
plt.legend(fontsize=9)


plt.figure(7,figsize=(5, 4))
plt.plot(xs,DELn,'+r',label='__nolegend__')
plt.plot(xs[-1],DELn[-1][0],'+r',label='Optimal damping - Linear')
plt.plot(xs,np.mean(DELn,1),'--r')
plt.plot([xs[0],xs[-1]],[1,1],':k')
hfont = {'fontname':'Arial'}
plt.xlabel('Mean (total) velocity at the reference height [m/s]',**hfont, fontsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.ylabel('Normalized DEL [-]',**hfont, fontsize=12)
plt.ylim([0.2,1.8])
plt.legend(fontsize=9)




