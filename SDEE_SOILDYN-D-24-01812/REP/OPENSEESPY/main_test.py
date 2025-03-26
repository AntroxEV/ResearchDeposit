# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:24:27 2023

@author: at924
"""
import pickle as pcle
import numpy as np
import scipy as sp
import OPS_15MW_GravityWind as ops

case ='_5ms'
target=188970.763000000

res=ops.modelWind(0.11436558952909699,target,case)


print(' value=' , res)


