# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:30:49 2024

@author: at924
"""

import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind_Random as ops

caseW ='_5ms'
ntask=1
QA='/MYmud' + caseW +'_SET' + str(ntask)+'.txt';
MYmud=np.loadtxt('./' + caseW + QA,skiprows=12001, usecols=1)
target=np.max(np.abs(MYmud))
print('target is', target)
ops.modelWind(0.25,target,caseW,ntask)