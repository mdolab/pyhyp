from pyhyp import pyHyp
fileName = 'm6_small.fmt'
#fileName = 'm6_large.fmt'
#fileName = 'm6_tiny.fmt'

options= {
    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile':fileName,
    
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65,
    's0':1.5e-5,
    'rMin':25,
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
    'writeMetrics': False,
    }

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('m6.cgns')
hyp.writePlot3D('m6.p3d')

