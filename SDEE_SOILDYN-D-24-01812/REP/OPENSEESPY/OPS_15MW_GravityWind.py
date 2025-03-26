# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 21:07:48 2023

@author: at924
"""
def modelWind(xdamp,target,tmax,caseW):

    print("=========================================================")
    
    print("IEA 15MW Reference Wind Turbine - Monopile - 2D - Gravity and Wind loads")
    
    import numpy as np
    #import matplotlib.pyplot as plt
    import openseespy.opensees as op
    import time
    
    startRun = time.time()
    
    #------------------------------------------------
    #TOOLS
    def CircularSection(D,t):
        A=np.pi*D**2/4-np.pi*(D-2*t)**2/4
        I=np.pi*D**4/64-np.pi*(D-2*t)**4/64
        Av=A*2/np.pi
        return A,I,Av
        
     #INITIAL DATA
     #INITIAL DATA
    if hasattr(xdamp,"__len__"):
        eta=xdamp[0]
    else:
        eta=xdamp
    print('eta',eta)
    #target=args[0]
    #caseW=args[1]
    #print('target =', target)
    #print('case =', caseW)
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
    
    #INPUT FILES
    
    FX=np.genfromtxt('./AerodynForces/FX'  + caseW + '.txt')
    FZ=np.genfromtxt('./AerodynForces/FZ' + caseW + '.txt')
    MY=np.genfromtxt('./AerodynForces/MY' + caseW + '.txt')
    #timeF=list(FX[:,0])
    windFX=list(FX[:,1])
    windFZ=list(FZ[:,1])
    windMY=list(MY[:,1])
    
    # SET UP ----------------------------------------------------------------------------
    op.wipe()						       # clear opensees model
    op.model('basic', '-ndm', 2, '-ndf', 3)	       # 2 dimensions, 3 dof per node
    
    # define GEOMETRY -------------------------------------------------------------
    # TOWER
    # nodal coordinates:
    op.node(1, 0., 15.)					   # node#, X Y
    op.node(2, 0., 28.)
    op.node(3, 0., 41.)
    op.node(4, 0., 54.)
    op.node(5, 0., 67.)					   
    op.node(6, 0., 80.)
    op.node(7, 0., 93.)
    op.node(8, 0., 106.)
    op.node(9, 0., 119.)					   
    op.node(10, 0., 132.)
    op.node(11, 0., 144.582,'-mass',28.28,28.28,0.) #TOP
    # TP
    # nodal coordinates:
    op.node(12, 0., 0.,'-mass', 100.00,100.00,125.00)					   
    op.node(13, 0., 5.)
    op.node(14, 0., 10)
    # FREE MONOPILE 
    # nodal coordinates:
    op.node(15, 0.,-30.)            #SB
    op.node(16, 0.,-25.)
    op.node(17, 0.,-20.)
    op.node(18, 0.,-15.)
    op.node(19, 0.,-10.)
    op.node(20, 0.,-5.)
    # RNA
    # nodal coordinates:
    op.node(21,-12.0515, 149.9858,'-mass',341.825,341.825,0.)           #HUB
    op.node(22,4.72,148.857,'-mass',646.895,646.895,0.)          #NACELLE
    # EMBEDDED MONOPILE 
    # nodal coordinates:
    for i in range(0,90):
        op.node(23+i, 0.,-75.+i/2)           #TOE
        #print('node',22+i,'z',-75+i/2)
    
    # define SOIL NODES -------------------------------------------------------
    # LEFT SIDE START FROM NODE 113
    # nodal coordinates:
    for i in range(0,90):
        op.node(113+i, -1.,-75.+i/2)           #TOE
        #print('LEFT node',112+i,'z',-75+i/2)
    
    # RIGHT SIDE START FROM NODE 203
    # nodal coordinates:
    for i in range(0,90):
        op.node(203+i, 1.,-75.+i/2)           #TOE
        #print('RIGHT node',202+i,'z',-75+i/2)
        
    # Boundary Conditions----------------------------------------------------
    # LEFT SIDE START FROM NODE 113
    for i in range(0,90):
        op.fix(113+i, 1, 1, 1) 			           # node DX DY RZ
    
    # RIGHT SIDE START FROM NODE 203
    for i in range(0,90):
        op.fix(203+i, 1, 1, 1) 			           # node DX DY RZ
    
    op.fix(23,0,1,0)                   #TOE ROLLER  
    
    # Constraints----------------------------------------------------
    op.rigidLink('beam', 11, 21)
    op.rigidLink('beam', 11, 22)
    
    # Define ELEMENTS -------------------------------------------------------------
    # define geometric transformation
    op.geomTransf('Corotational', 1)  		       # associate a tag to transformation
    #geomTransf('Linear', 1)  		       # associate a tag to transformation
    
    # Define CONNECTIVITY --------------------------------------------------------
    # element('elasticBeamColumn', eleTag, *eleNodes, Area, E_mod, Iz, transfTag, <'-mass', mass>, <'-cMass'>, <'-release', releaseCode>)
     #element('ElasticTimoshenkoBeam', eleTag, *eleNodes, E_mod, G_mod, Area, Iz, Avy, transfTag, <'-mass', massDens>, <'-cMass'>)
    # TOWER
    for i in range(0,10):
        A,I,Av=CircularSection(Dtw[i],ttw[i])
        #element('ElasticTimoshenkoBeam', 1+i, 1+i, 2+i, E, G,A,I, Av,1,'-mass',A*rho)
        op.element('elasticBeamColumn', 1+i, 1+i, 2+i, A, E,I,1,'-mass',A*rho)
    
    # TP
    for i in range(0,2):
        A,I,Av=CircularSection(Dtp[i],ttp[i])
        #element('ElasticTimoshenkoBeam', 12+i, 12+i, 13+i, E, G,A,I, Av,1,'-mass',A*rho)
        op.element('elasticBeamColumn', 12+i, 12+i, 13+i, A, E,I,1,'-mass',A*rho)
    
    A,I,Av=CircularSection(Dtp[2],ttp[2])
    #element('ElasticTimoshenkoBeam', 12+2, 14, 1, E, G, A,I, Av,1,'-mass',A*rho,)
    op.element('elasticBeamColumn', 12+2, 14, 1, A, E, I,1,'-mass',A*rho,)
    # FREE MONOPILE 
    for i in range(0,5):
        A,I,Av=CircularSection(Dfmp[i],tfmp[i])
        #element('ElasticTimoshenkoBeam', 15+i, 15+i, 16+i, E, G,A, I, Av,1,'-mass',A*rho)
        op.element('elasticBeamColumn', 15+i, 15+i, 16+i, A, E, I,1,'-mass',A*rho)
    	
    A,I,Av=CircularSection(Dfmp[5],tfmp[5])
    #element('ElasticTimoshenkoBeam', 15+5, 15+5, 12, E, G,A, I, Av,1,'-mass',A*rho) 
    op.element('elasticBeamColumn', 15+5, 15+5, 12, A, E, I,1,'-mass',A*rho)
    # EMBEDDED MONOPILE 
    A,I,Av=CircularSection(Demp,temp)
    for i in range(0,89):
        #element('ElasticTimoshenkoBeam', 21+i, 23+i, 24+i, E, G,A, I, Av,1,'-mass',A*rho)
        op.element('elasticBeamColumn', 21+i, 23+i, 24+i, A, E, I,1,'-mass',A*rho)
    #element('ElasticTimoshenkoBeam', 21+89, 23+89, 15, E, G,A, I, Av,1,'-mass',A*rho)
    op.element('elasticBeamColumn', 21+89, 23+89, 15, A, E, I,1,'-mass',A*rho)	
    # SPRINGS 
    Es=[37659.06417,38595.36532,39492.07176,40349.18346,41166.70045,41944.62271,42682.95024,43381.68306,44040.82114,44660.36451,45240.31315,45780.66706,46281.42625,46742.59072,47164.16047,47546.13549,47888.51578,48191.30135,48454.4922,	48678.08832,48862.08972,49006.4964,	49111.30835,49176.52557,49202.14808,49188.17585,49134.60891,49041.44724,48908.69084,48736.33973,48524.39389,48272.85332,47981.71803,47650.98802,47280.66328,46870.743812,46421.22963,45932.12072,45403.41708,44835.11872,44227.22564,43579.73783,42892.6553,42165.97805,41399.70607,40593.83937,39748.37794,38863.32179,37938.67091,36974.42531,35970.58499,34927.14994,34613.30472,35058.74537,35504.18601,	35949.62666,	36395.0673,	36840.50795,	37285.9486,	37731.38924,	38176.82989,	38622.27053,	39067.71118,	39513.15183,	39958.59247,	40404.03312,	40849.47376,	41294.91441,	41740.35505,	272064.0605,	274936.7902,	277809.52,280682.2498,283554.9795,286427.7093,289300.439,	292173.1688,295045.8986,297918.6283,300791.3581,303664.0878,306536.8176,309409.5474,312282.2771,315155.0069,318027.7367,320900.4664,323773.1962,326645.9259,303177.084]
    Ers=[2114250.293206,	2166815.991226,	2217158.767384,	2265278.621681,	2311175.554117,	2354849.564692,	2396300.653407,	2435528.820260,	2472534.065252,	2507316.388383,	2539875.789653,	2570212.269062,	2598325.826610,	2624216.462297,	2647884.176123,	2669328.968088,	2688550.838191,	2705549.786434,	2720325.812816,	2732878.917337,	2743209.099997,	2751316.360795,	2757200.699733,	2760862.116810,	2762300.612025,	2761516.185380,	2758508.836874,	2753278.566506,	2745825.374278,	2736149.260188,	2724250.224238,	2710128.266426,	2693783.386754,	2675215.585220,	2654424.861825,	2631411.216570,	2606174.649453,	2578715.160475,	2549032.749637,	2517127.416937,	7725583.905764,	7612481.143679,	7492462.000141,	7365526.475151,	7231674.568708,	7090906.280813,	6943221.611465,	6788620.560665,	6627103.128412,	6458669.314707,	6283319.119549,	6101052.542938,	6046230.257259,	6124039.548600,	6201848.839941,	6279658.131282,	6357467.422623,	6435276.713964,	6513086.005305,	6590895.296646,	6668704.587987,	6746513.879328,	6824323.170669,	6902132.462010,	6979941.753351,	7057751.044692,	7135560.336033,	7213369.627374,	7291178.918715,	7368988.210056,	7446797.501397,	7524606.792738,	7602416.084079,	7680225.375420,	7758034.666761,	7835843.958102,	7913653.249443,	7991462.540784,	8069271.832125,	8147081.123466,	8224890.414807,	8302699.706148,	8380508.997489,	8458318.288830,	8536127.580171,	8613936.871511,	8691746.162852,	8769555.454193,	8847364.745534,9.9854E+006]
    # element('twoNodeLink', eleTag, *eleNodes, '-mat', *matTags, '-dir', *dir, <'-orient', *vecx, *vecyp>, <'-pDelta', *pDeltaVals>, <'-shearDist', *shearDist>, <'-doRayleigh'>, <'-mass', m>)
    for i in range(0,90):
        op.uniaxialMaterial('Elastic', 1+i, Es[i], np.multiply(eta,Es[i]))
        op.uniaxialMaterial('Elastic', 91+i, Ers[i], np.multiply(eta,Ers[i]))
    
    # LEFT SIDE START FROM NODE 113
    # nodal coordinates:
    for i in range(0,90):
        op.element('twoNodeLink', 111+i, 113+i, 23+i, '-mat', 90-i, 180-i, '-dir', 1,3)
    # RIGHT SIDE START FROM NODE 203
    for i in range(0,90):
        op.element('twoNodeLink', 201+i, 203+i, 23+i, '-mat', 90-i, 180-i, '-dir', 1,3)
    
    # Define ANALYSIS -------------------------------------------------------------
    
    
    # define GRAVITY ANALYSIS------------------------------------------------------
    
    op.timeSeries('Constant', 1, '-factor', 9.81) 
    op.pattern('UniformExcitation', 1, 2,'-accel', 1)
    op.constraints('Transformation')  				# how it handles boundary conditions
    op.numberer('Plain')			    # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral')		    # how to store and solve the system of equations in the analysis
    op.algorithm('Newton')                 # use Linear algorithm for linear analysis
    op.integrator('LoadControl', 0.1)    # determine the next time step for an analysis
    op.analysis('Static')   # define type of analysis: time-dependent
    uTOP = [0.0]
    MYmud=[0.0]
    op.analyze(1,0.1)
    uTOP.append(op.nodeDisp(11,1))
    MYmud.append(op.eleForce(15, 3))
    print('GRAVITY: uTOP = ', uTOP, '\n')
    op.loadConst('-time', 0.0)				# hold gravity constant and restart time
    # define WIND ANALYSIS------------------------------------------------------
    op.timeSeries('Path', 2,   '-dt', 0.005, '-values', *windFX,'-prependZero', '-factor',0.001)
    op.timeSeries('Path', 3,  '-dt', 0.005, '-values', *windFZ,'-prependZero', '-factor',0.001)
    op.timeSeries('Path', 4,  '-dt', 0.005, '-values', *windMY,'-prependZero', '-factor',0.001)
     
    op.wipeAnalysis()
    
    #region(1,'-ele', range(1,111), '-rayleigh', 0.0348, 0., 0.0054, 0.)
    op.rayleigh(0.0348, 0.0054, 0., 0.) #0.02  damping
    
    op.pattern('Plain', 2, 2)
    op.load(21, 1,0,0)
    op.pattern('Plain', 3, 3)
    op.load(21, 0,1,0)
    op.pattern('Plain', 4, 4)
    op.load(21, 0,0,-1)
    
    
    op.constraints('Transformation')  				# how it handles boundary conditions
    op.numberer('Plain')			    # renumber dof's to minimize band-width (optimization), if you want to
    op.system('BandGeneral')		    # how to store and solve the system of equations in the analysis
    op.test('NormDispIncr', 1.0e-4,  3 )
    op.algorithm('Newton')                 # use Linear algorithm for linear analysis
    op.integrator('HHT', 0.9)    # determine the next time step for an analysis
    op.analysis('Transient')   # define type of analysis: time-dependent
    # SOLVER GRAVITY -------------------------------------------------------------
    tCurrent = op.getTime()
    
    Atime = [0.0, tCurrent]
    aTOP = [0.0, 0.0]
    
    while tCurrent <= tmax:
        ok = op.analyze(1, .01)
        if ok != 0:
            print('noncovergence')
        tCurrent = op.getTime()
        Atime.append(tCurrent)
        uTOP.append(op.nodeDisp(11,1))
        aTOP.append(op.nodeAccel(11,1))
        MYmud.append(op.eleForce(15, 3))
    
        
    
    #plt.plot(uTOP)
    #print('time',  time, "uTOP = \n", uTOP)
    #loadConst('-time', 0.0)				# hold gravity constant and restart time

    #np.savetxt('uTOP' + caseW + '_OPT.txt',uTOP)
    #np.savetxt('aTOP'+ caseW  + '_OPT.txt',aTOP)
    #np.savetxt('time'+ caseW + '_OPT.txt',Atime)
    #np.savetxt('MYmud'+ caseW + '_OPT.txt',MYmud)
    np.savetxt('uTOP' + caseW + '_OPTVER.txt',uTOP)
    np.savetxt('aTOP'+ caseW  + '_OPTVER.txt',aTOP)
    np.savetxt('time'+ caseW + '_OPTVER.txt',Atime)
    np.savetxt('MYmud'+ caseW + '_OPTVER.txt',MYmud)
    #np.savetxt('uTOP' + caseW + '.txt',uTOP)
    #np.savetxt('aTOP'+ caseW  + '.txt',aTOP)
    #np.savetxt('time'+ caseW + '.txt',Atime)
    #np.savetxt('MYmud'+ caseW + '.txt',MYmud)
    optim=np.divide(np.abs(np.max(np.abs(MYmud[12000:]))-target),target)
    print('optim',optim)
    print('max MYmud',np.max(np.abs(MYmud[12000:])))
    elapsed_time = time.time() - startRun
    print('elapsed time (minutes):',elapsed_time/60)
    
    return optim
