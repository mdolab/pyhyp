from pyhyp import pyHyp

options= {
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65,
    's0':1.5e-6,
    'rMin':60,

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
    'epsI': 2.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.005,
    'volSmoothIter': 10,
    
    # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-5,
    'kspMaxIts': 500,
    'preConLag': 5,
    'kspSubspaceSize':50,
    }

hyp = pyHyp.pyHyp('3d',fileName='m6_small.fmt', options=options, zMirror=True)
#hyp = pyHyp.pyHyp('3d',fileName='m6_large.fmt', options=options, zMirror=True)

hyp.run()
hyp.writeCGNS('m6.cgns')

