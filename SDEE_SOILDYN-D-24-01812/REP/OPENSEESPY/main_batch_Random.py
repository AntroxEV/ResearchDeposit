# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""

import OPS_15MW_GravityWind_Random as ops
from multiprocessing import Pool

NUMBER_OF_TASKS = 2
caseW ='_25ms'

target= 315386.935000000    #5ms: 188970.763000000 (DEGR:163202.628000000) 10ms:521309.100000000  12_5ms: 398961.227000000 15ms: 360037.651000000 15.5ms: 345150.785000000 20ms: 338790.933000000  (DEGR:292670.165000000 ) 25ms: 315386.935000000 50ms:404947.587000000 (250041.192000000) 
eta=[[9.301757812499847344e-02],[8.837402343750002665e-01]]
anal=[1,3]
def RUN(ntask): 
    ops.modelWind(eta[ntask],target, 720,caseW,anal[ntask] )
    print('DONE eta',eta[ntask],' case',anal[ntask])


if __name__ == '__main__':
    #fpath="./15ms"
    pool = Pool(NUMBER_OF_TASKS)
    ntask = [x for x in range(0, NUMBER_OF_TASKS)]
    result = pool.map(RUN,ntask)
    print(result)
    pool.close()
    pool.join()
    
