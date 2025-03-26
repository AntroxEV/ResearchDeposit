# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind_DEL_Random as ops
import rainflow as rf
from multiprocessing import Pool

NUMBER_OF_TASKS = 10
SETS=10
caseW ='_25ms'
guessd=0.2

def STOP(intermediate_result):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.01:
        raise StopIteration

    

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
    optis={'maxiter':15,'disp':True,'xatol':1e-4,'return_all':True,'fatol':0.1}
    res = sp.minimize(ops.modelWind, guessd,method='Nelder-Mead',bounds=sp.Bounds(0., 10,False),callback=STOP,args=(target,tmax,caseW,ntask),options=optis)
    file_pi = open('OptRes'+caseW+'_SET' + str(ntask)+'_DEL.obj', 'wb')
    pcle.dump(res, file_pi)
    tmp=res.x
    print('optimal value=' , tmp[0])
    file_pi.close()
    np.savetxt('OptRes'+caseW+'_SET' + str(ntask)+'_DEL.txt', (tmp[0],res.fun))
    return 0



def RUN(ntask): 
    Nref=2E8;
    mcoef=4;
    rngT=[]
    countT=[]
    QA='/MYmud' + caseW +'_SET' + str(ntask)+'.txt';
    MYmud=np.loadtxt('./' + caseW + QA,skiprows=12001, usecols=1)
    #array_out = rf.rainflow(MYmud)
    #array_out = array_out[:,array_out[0,:].argsort()]
    #target=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
    for rng, a , count, b,c in rf.extract_cycles(MYmud):
        rngT.append(rng)
        countT.append(count)
    target = np.power(np.divide(np.sum(np.multiply(np.power(rngT,mcoef),countT)),Nref),1/mcoef)
    print('target DEL',target)
    OPTIM1b(target,720,caseW,ntask)
    
    
if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(1, SETS+1)]
    #ntask = [5,6,7,8,9]
    #RUN(ntask[0])
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()