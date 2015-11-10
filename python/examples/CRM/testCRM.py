from pyhyp import pyHyp

options= {

    # ---------------------------
    #        Input File
    # ---------------------------
    'inputFile':'dpw.fmt',

    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    'N': 65, 
    's0':1e-3,
    'gridRatio':1.20,
    'rMin':5,
    'mirror':'y',
    'symTol':1e-6,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    'ps0':1e-3,
    'pGridRatio':1.15,

    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    'epsE': 2.0,
    'epsI': 4.0,
    'theta': 3.0,
    'volCoef': .16,
    'volBlend': 0.001,
    'volSmoothIter': 10,

   # ---------------------------
    #   Solution Parameters
    # ---------------------------
    'kspRelTol': 1e-10,
    'kspMaxIts': 500,
    'preConLag': 10,
    'kspSubspaceSize':50,
    }

hyp = pyHyp(options=options, debug=True)
hyp.run()
hyp.writeCGNS('dpw.cgns')

