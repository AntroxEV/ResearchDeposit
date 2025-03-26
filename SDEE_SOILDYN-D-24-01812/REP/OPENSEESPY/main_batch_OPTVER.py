# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import OPS_15MW_GravityWind as ops
from multiprocessing import Pool
import numpy as np


NUMBER_OF_TASKS = 5
caseW =['_5ms','_10ms','_12_5ms','_15ms','_17_5ms','_20ms','_25ms','_30ms','_40ms','_50ms']
xs=[5.,10.,12.5,15.,17.5,20.,25,30,40,50]
res3=np.array([ 0.10527403, -0.01943716,  0.00194243]) #OPT 
slope=0.027638125000000006
intercept=-0.03637812500000004

print('total number of analyses',len(caseW))

def RUN(ntask):
    case='./SEISMO_NL/MYmud' + caseW[ntask] + '.txt'
    AA=np.genfromtxt(case)
    #eta=slope*xs[ntask]+intercept
    target=np.max(np.abs(AA[12000:,1]))
    print('target MYmud',target) 
    eta=res3[2]*xs[ntask]**2+res3[1]*xs[ntask]+res3[0]
    optim=ops.modelWind(eta,target,720,caseW[ntask]  )
    print('DONE',caseW[ntask],'optim',optim)


if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(0, len(caseW))]
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()
    
