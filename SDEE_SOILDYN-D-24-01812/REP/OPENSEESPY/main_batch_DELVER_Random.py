# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import OPS_15MW_GravityWind_DEL_Random as ops
from multiprocessing import Pool
import numpy as np
import itertools as it


NUMBER_OF_TASKS = 20
SETS=range(1,11)
caseW =np.array(['_5ms','_10ms','_15ms','_20ms','_25ms','_50ms'])


xs=[5.,10.,15.,20.,25,50]
res3=np.array([-0.22822007,  0.05999973, -0.00107024])#DEL

print('total number of analyses',len(caseW))

def RUN(param):
    ntask=param[1]
    casew=param[0]
    index=np.where(caseW == casew)
    vel=xs[index[0][0]]
    eta=res3[2]*vel**2+res3[1]*vel+res3[0]
    ops.modelWind(eta,1e6,720,casew,ntask  )
    print('DONE',casew,'SET',ntask)


if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    paramlist = list(it.product(caseW,SETS))
    result = pool.map(RUN,paramlist)
    print(result)
    pool.close()
    pool.join()
    
