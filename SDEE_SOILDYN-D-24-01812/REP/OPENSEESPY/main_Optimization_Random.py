# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind_Random as ops
import subprocess
from multiprocessing import Pool

NUMBER_OF_TASKS = 10
caseW ='_50ms'

def STOP(intermediate_result: sp.OptimizeResult):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.01:
        raise StopIteration

    
    

def OPTIM1b(target,tmax,caseW,ntask):
    optis={'maxiter':15,'disp':True,'xatol':1e-4,'return_all':True,'fatol':0.1}
    res = sp.minimize(ops.modelWind, 3.9,method='Nelder-Mead',bounds=sp.Bounds(0., 10,False),callback=STOP,args=(target,tmax,caseW,ntask),options=optis)
    file_pi = open('OptRes'+caseW+'_SET' + str(ntask)+'_OPT.obj', 'wb')
    pcle.dump(res, file_pi)
    tmp=res.x
    print('optimal value=' , tmp[0])
    file_pi.close()
    np.savetxt('OptRes'+caseW+'_SET' + str(ntask)+'.txt', (tmp[0],res.fun))
    return 0

def RUN(ntask): 
    QA='/MYmud' + caseW +'_SET' + str(ntask)+'.txt';
    MYmud=np.loadtxt('./' + caseW + QA,skiprows=12001, usecols=1)
    target=np.max(np.abs(MYmud))
    print('target is', target, 'at SET',str(ntask))
    OPTIM1b(target,720,caseW,ntask)
    
    
if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(1, NUMBER_OF_TASKS+1)]
    #RUN(ntask[0])
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()