from pyhypellip import pyHyp

options= {
    # ---------------------------
    #        Input Information
    # ---------------------------
    'dimension':'3d',
    #'inputFile':'m6_small.fmt',
    #'inputFile':'m6_large.fmt',
    'inputFile':'m6_tiny.fmt',
    'mode':'hyperbolic',
    
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 33,
    's0':1.5e-4,
    'rMin':1,
    'mirror':'z',
    
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':1.5e-6,
    'pGridRatio':1.1,
    'cMax':6,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 4.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.005,
    'volSmoothIter': 100,
    
    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 500,
    'preConLag': 10,
    'kspSubspaceSize':50,
    'writeMetrics': True,
    }

hyp = pyHyp(options=options)
#hyp.writeFEMesh('test.dat')
hyp.run()
# hyp.hyp.run3delliptic()
hyp.writeCGNS('m6.cgns')
hyp.writePlot3D('m6.p3d')

