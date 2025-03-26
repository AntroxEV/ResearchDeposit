# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind_Random as ops
import OPS2_15MW_GravityWind_Random as ops2
import subprocess
from multiprocessing import Pool


caseW ='_5ms'

def STOP(intermediate_result: sp.OptimizeResult):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.005:
        return StopIteration
    else:
        return 0
    
    
    

def OPTIM1(target,caseW,ntask):
    optis={'maxiter':15,'disp':3,'xatol':1e-3}
    res = sp.minimize_scalar(ops.modelWind, bounds=(0.25, 0.9),args=(target,caseW,ntask), method='bounded',options=optis)
    file_pi = open('OptRes'+caseW+'_SET' + str(ntask)+'.obj', 'wb')
    pcle.dump(res, file_pi)
    print('optimal value=' , res.x)
    file_pi.close()
    np.savetxt('OptRes'+caseW+'_SET' + str(ntask)+'.txt', (res.x,res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0
    
def OPTIM1b(target,tmax,caseW,ntask):
    optis={'maxiter':15,'disp':True,'xatol':1e-4,'return_all':True,'fatol':5e-3}
    res = sp.minimize(ops.modelWind, 0.75,method='Nelder-Mead',bounds=sp.Bounds(0., 0.99,False),args=(target,tmax,caseW,ntask),options=optis)
    file_pi = open('OptRes'+caseW+'_SET' + str(ntask)+'.obj', 'wb')
    pcle.dump(res, file_pi)
    print('optimal value=' , res.x)
    print('start final analysis')
    ops.modelWind(res.x,target,tmax,caseW,ntask)
    file_pi.close()
    np.savetxt('OptRes'+caseW+'_SET' + str(ntask)+'.txt', (res.x,res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0

def OPTIM2(target,caseW):
    optis={'maxiter':15,'disp':3,'xatol':1e-4}
    res = sp.minimize_scalar(ops2.modelWind, bounds=(0.0, 0.9),args=(target,caseW,ntask), method='bounded',options=optis)
    file_pi = open('OptRes'+caseW+'_SET' + str(ntask)+'.obj', 'wb')
    pcle.dump(res, file_pi)
    print('optimal value=' , res.x)
    file_pi.close()
    np.savetxt('OptRes'+caseW+'_SET' + str(ntask)+'.txt', (res.x,res.fun))
    return 0

def RUN(ntask): 
    QA='/MYmud' + caseW +'_SET' + str(ntask)+'.txt';
    MYmud=np.loadtxt('./' + caseW + QA,skiprows=12001, usecols=1)
    target=np.max(np.abs(MYmud))
    indexmax=np.argmax((np.abs(MYmud)))
    print('max at:',indexmax, 'at SET',str(ntask))
    print('target is', target, 'at SET',str(ntask))
    tmax=(indexmax+12001+10)*0.01
    OPTIM1b(target,tmax,caseW,ntask)
    
    
if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(2)
    ntask = [4,6]
    #RUN(ntask[0])
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()