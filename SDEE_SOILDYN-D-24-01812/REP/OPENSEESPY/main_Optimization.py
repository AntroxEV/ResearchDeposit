# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy.optimize as sp
import OPS_15MW_GravityWind as ops
import OPS2_15MW_GravityWind as ops2

caseW ='_12_5ms'
target= 398961.227000000     #5ms: 188970.763000000 (DEGR:163202.628000000) 10ms:521309.100000000  12_5ms: 398961.227000000 15ms: 360037.651000000 15.5ms: 345150.785000000 20ms: 338790.933000000  (DEGR:292670.165000000 ) 25ms: 315386.935000000 50ms:404947.587000000 (250041.192000000) 

def STOP(intermediate_result):
    print('intermediate_results',intermediate_result)
    if intermediate_result.fun < 0.005:
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
    res = sp.minimize(ops.modelWind, 1.2,method='Nelder-Mead',bounds=sp.Bounds(0., 10.,False),callback=STOP,args=(target,720,caseW),options=optis)
    file_pi = open('OptRes'+caseW+'.obj', 'wb')
    pcle.dump(res, file_pi)
    tmp=res.x
    print('optimal value=' , tmp[0])
    file_pi.close()
    #print('start final analysis')
    #ops.modelWind(res.x,target,720,caseW)
    np.savetxt('OptRes'+caseW+'.txt', (tmp[0],res.fun))
    #filehandler = open(filename, 'r') 
    #object = pickle.load(filehandler)
    return 0


if __name__ == '__main__':
    OPTIM1b(target,caseW)
    
