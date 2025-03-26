# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind_DEL as ops
import OPS2_15MW_GravityWind as ops2
import rainflow as rf
from multiprocessing import Pool

NUMBER_OF_TASKS = 3





def STOP(intermediate_result):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.01:
        raise StopIteration


def OPTIM(ntask):
    Nref=2E8;
    mcoef=4;
    caseW =['_5ms','_10ms','_15ms','_20ms','_25ms']
    x0=[0,0.2,0.78,0.2,0.65]
    case='./SEISMO_NL/MYmud' + caseW[ntask] + '.txt'
    AA=np.genfromtxt(case)
    array_out = rf.rainflow(AA[12000:,1])
    array_out = array_out[:,array_out[0,:].argsort()]
    target=np.power(np.divide(np.sum(np.multiply(np.power(array_out[0],mcoef),array_out[3])),Nref),1/mcoef)
    print('target DEL',target)
    optis={'maxiter':15,'disp':True,'xatol':1e-4,'return_all':True,'fatol':5e-3}
    res = sp.minimize(ops.modelWind, x0[ntask],method='Nelder-Mead',bounds=sp.Bounds(0., 10.,False),callback=STOP,args=(target,720,caseW[ntask]),options=optis)
    file_pi = open('OptRes'+caseW+'_DEL.obj', 'wb')
    pcle.dump(res, file_pi)
    tmp=res.x
    print('optimal value=' , tmp[0])
    file_pi.close()
    #print('start final analysis')
    #ops.modelWind(res.x,target,720,caseW)
    np.savetxt('OptRes'+caseW+'_DEL.txt', (tmp[0],res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0


if __name__ == '__main__':
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(0, 5)]
    result = pool.map(OPTIM,ntask)
    print(result)
    pool.close()
    pool.join()

    
