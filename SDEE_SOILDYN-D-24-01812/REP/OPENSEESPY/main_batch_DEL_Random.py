# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import OPS_15MW_GravityWind_DELsens as ops
from multiprocessing import Pool
import numpy as np
import rainflow as rf

NUMBER_OF_TASKS = 6
caseW ='_50ms'

Nref=2E8;
mcoef=4;
case='./SEISMO_NL/MYmud' + caseW + '.txt'
AA=np.genfromtxt(case)
array_out = rf.rainflow(AA[12000:,1])
array_out = array_out[:,array_out[0,:].argsort()]
target=np.power(np.divide(np.sum(np.power(np.multiply(array_out[0],array_out[3]),mcoef)),Nref),1/mcoef)
print('target DEL',target)

eta=[0.,0.2,0.4,0.6,0.8,3.]

def RUN(ntask): 
    penal=ops.modelWind(eta[ntask],target,720,caseW )
    print('DONE eta',eta[ntask],'penal',penal)


if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(0, NUMBER_OF_TASKS)]
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()
    
