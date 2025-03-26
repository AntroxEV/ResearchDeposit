# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 21:07:48 2023

@author: at924
"""

print("=========================================================")

print("IEA 15MW Reference Wind Turbine - Monopile - 2D - Gravity and Wind loads")

import numpy as np
from openseespy.opensees import *
#TOOLS
def CircularSection(D,t):
    A=np.pi*D**2/4-np.pi*(D-2*t)**2/4
    I=np.pi*D**4/64-np.pi*(D-2*t)**4/64
    Av=A*2/np.pi
    return A,I,Av
    
 #INITIAL DATA
E=2E8
v=0.3
rho=8.6
G=E/2/(1+v)
Dtw=[10.,10.,9.926,9.443,8.833,8.151,7.39,6.909,6.748,6.572] 
ttw=[0.03950,0.03646,0.03646,0.03378,0.03219,0.03071,0.02721,0.02083,0.02083,0.02400,]
Dtp=[10,10,10]
ttp=[0.04326,0.04324,0.04106]
Dfmp=[10]*6
tfmp=[0.05534,0.05345,0.05151,0.04953,0.04752,0.04552]
Demp=10
temp=0.05534




# SET UP ----------------------------------------------------------------------------
wipe()						       # clear opensees model
model('basic', '-ndm', 2, '-ndf', 3)	       # 2 dimensions, 3 dof per node

# define GEOMETRY -------------------------------------------------------------
# TOWER
# nodal coordinates:
node(1, 0., 15.)					   # node#, X Y
node(2, 0., 28.)
node(3, 0., 41.)
node(4, 0., 54.)
node(5, 0., 67.)					   
node(6, 0., 80.)
node(7, 0., 93.)
node(8, 0., 106.)
node(9, 0., 119.)					   
node(10, 0., 132.)
node(11, 0., 144.582,'-mass',28.28,28.28,0.) #TOP
# TP
# nodal coordinates:
node(12, 0., 0.,'-mass', 100.00,100.00,125.00)					   
node(13, 0., 5.)
node(14, 0., 10)
# FREE MONOPILE 
# nodal coordinates:
node(15, 0.,-30.)            #SB
node(16, 0.,-25.)
node(17, 0.,-20.)
node(18, 0.,-15.)
node(19, 0.,-10.)
node(20, 0.,-5.)
# RNA
# nodal coordinates:
node(21,-12.0515, 149.9858,'-mass',341.825,341.825,0.)           #HUB
node(22,4.72,148.857,'-mass',646.895,646.895,0.)          #NACELLE

    
# Boundary Conditions----------------------------------------------------

fix(15,1,1,1)                   #TOE ROLLER  



# Constraints----------------------------------------------------
rigidLink('beam', 11, 21)
rigidLink('beam', 11, 22)

# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation
geomTransf('Linear', 1)  		       # associate a tag to transformation
# geomTransf('PDelta', 1)  		       # associate a tag to transformation

# connectivity:
# element('elasticBeamColumn', eleTag, *eleNodes, Area, E_mod, Iz, transfTag, <'-mass', mass>, <'-cMass'>, <'-release', releaseCode>)
 #element('ElasticTimoshenkoBeam', eleTag, *eleNodes, E_mod, G_mod, Area, Iz, Avy, transfTag, <'-mass', massDens>, <'-cMass'>)
# TOWER
for i in range(0,10):
    A,I,Av=CircularSection(Dtw[i],ttw[i])
    element('ElasticTimoshenkoBeam', 1+i, 1+i, 2+i, E, G,A,I, Av,1,'-mass',A*rho)	
# TP
for i in range(0,2):
    A,I,Av=CircularSection(Dtp[i],ttp[i])
    element('ElasticTimoshenkoBeam', 12+i, 12+i, 13+i, E, G,A,I, Av,1,'-mass',A*rho)
A,I,Av=CircularSection(Dtp[2],ttp[2])
element('ElasticTimoshenkoBeam', 12+2, 14, 1, E, G, A,I, Av,1,'-mass',A*rho)
# FREE MONOPILE 
for i in range(0,5):
    A,I,Av=CircularSection(Dfmp[i],tfmp[i])
    element('ElasticTimoshenkoBeam', 15+i, 15+i, 16+i, E, G,A, I, Av,1,'-mass',A*rho)	
A,I,Av=CircularSection(Dfmp[5],tfmp[5])
element('ElasticTimoshenkoBeam', 15+5, 15+5, 12, E, G,A, I, Av,1,'-mass',A*rho) 

# Define ANALYSIS -------------------------------------------------------------

wipeAnalysis()			     # clear previously-define analysis parameters
constraints('Transformation')    	 # how it handles boundary conditions
#constraints('Lagrange', 1.0, 1.0)
numberer('Plain')    # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral') # how to store and solve the system of equations in the analysis
numEigen = 10
eigenValues = eigen(numEigen)
natfreq=[(eigenValues[i]**0.5)/2/np.pi for i in range(0,numEigen)]
print("eigen values:",natfreq)
  