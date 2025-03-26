#!/usr/bin/env python
"""
Demonstration of rainflow cycle counting

Contact: Jennifer Rinker, Duke University
Email:   jennifer.rinker-at-duke.edu
"""
import matplotlib.pyplot as plt
import rainflow as rf
import numpy as np




caseW='_10ms'

case='./SEISMO_NL/MYmud' + caseW + '.txt'
AA=np.genfromtxt(case)
plt.figure(1,figsize=(4, 4))
plt.plot(AA[12000:,0]-120,AA[12000:,1],'k',label='Nonlinear Model',linewidth=0.5)

array_out = rf.rainflow(AA[12000:,1])
# sort array_out by cycle range
array_out = array_out[:,array_out[0,:].argsort()]
    


# =============================================================================
# case='./SEISMO_L/DAMP0/MYmud' + caseW + '.txt'
# BB=np.genfromtxt(case) #5ms-10ms-15ms-20ms-25ms gravity
# plt.plot(BB[12000:,0]-120,BB[12000:,1],':r',label='Linear Model',linewidth=0.5)
# plt.legend()
# 
# plt.figure(2,figsize=(4, 4))
# plt.plot(AA[12000:,0]-120,AA[12000:,1],'k',label='Nonlinear Model',linewidth=0.5)
# case="C:/Users/at924/OneDrive - University of Exeter/Research/EPSRC Project/WP3 (WU)-Design Formulas and Simplified SSI model/Damping linearisation/Scripts_Optimization/Results/NODEGR/MYmud" + caseW + '_OPT.txt'
# CC=np.genfromtxt(case) #5ms-10ms-15ms-20ms-25ms gravity
# plt.plot(np.multiply(range(0,len(CC[12001:])),0.01),CC[12000:-1],':b',label='Optimized Model',linewidth=0.5)
# plt.legend()
# 
# 
# =============================================================================

