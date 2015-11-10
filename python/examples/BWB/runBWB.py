from pyhyp import pyHyp
fileName = 'bwb.fmt'
options= {
    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile':fileName,

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65,
    's0':4e-6,
    'rMin':25,
    'mirror':'z',
    
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':5e-6,
    'pGridRatio':1.1,
    'cMax':6,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 1.0,
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.0001,
    'volSmoothIter': 20,
    
    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    'writeMetrics': False,
    }

hyp = pyHyp(options=options)
hyp.run()
hyp.writeCGNS('bwb.cgns')

