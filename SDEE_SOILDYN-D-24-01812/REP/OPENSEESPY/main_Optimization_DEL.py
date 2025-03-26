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

caseW ='_50ms'

Nref=2E8;
mcoef=4;
case='./SEISMO_NL/MYmud' + caseW + '.txt'
AA=np.genfromtxt(case)
array_out = rf.rainflow(AA[12000:,1])
array_out = array_out[:,array_out[0,:].argsort()]
target=np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef)
print('target DEL',target)


def STOP(intermediate_result):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.01:
        raise StopIteration


def OPTIM1(target,caseW):
    optis={'maxiter':15,'disp':3,'xatol':1e-3}
    res = sp.minimize_scalar(ops.modelWind, bounds=(0.4, 1.999),args=(target,caseW), method='bounded',options=optis)
    file_pi = open('OptRes'+caseW+'.obj', 'wb')
    pcle.dump(res, file_pi)
    print('optimal value=' , res.x)
    file_pi.close()
    np.savetxt('OptRes'+caseW+'.txt', (res.x,res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0

def OPTIM2(target,caseW):
    optis={'maxiter':15,'disp':True,'xtol':1e-3,'ftol':0.1}
    res = sp.minimize(ops2.modelWind, [0.35,0.35],args=(target,caseW),bounds=[(0.2, 0.9),(0.2,0.9)], method='Powell',options=optis)
    file_pi = open('OptRes2'+caseW+'.obj', 'wb')
    pcle.dump(res, file_pi)
    print('optimal value=' , res.x)
    file_pi.close()
    #np.savetxt('OptRes2'+caseW+'.txt', (res.x,res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0


def OPTIM1b(target,caseW):
    optis={'maxiter':15,'disp':True,'xatol':1e-4,'return_all':True,'fatol':5e-3}
    res = sp.minimize(ops.modelWind, 0.1,method='Nelder-Mead',bounds=sp.Bounds(0., 10.,False),callback=STOP,args=(target,720,caseW),options=optis)
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
    OPTIM1b(target,caseW)
    
